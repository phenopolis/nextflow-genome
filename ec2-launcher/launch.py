import os
import argparse
import yaml
import time
import hashlib
import logging
from datetime import datetime
import boto3

COMPONENTS = {
    1: 'align',
    2: 'variantCall',
    3: 'annotation',
}


def get_report_s3_path(params):
    now = datetime.today().strftime('%Y-%m-%d-%H-%M-%S')
    randomiser = hashlib.sha224(params['input_table'].encode('utf-8')).hexdigest()
    return f"{params['report_path']}/{now}-{randomiser}"

def main(params):
    # configure logger
    os.makedirs('./log', exist_ok=True)
    now = datetime.today().strftime('%Y-%m-%d-%H-%M-%S')
    task_name = f"{now}-{params['input_table'].split('/')[-1]}"
    logfile = os.path.join('./log', f"{task_name}.log") 
    logging.basicConfig(filename=logfile, level=logging.DEBUG)
    print(f"log file: {logfile}")

    # report s3 path
    report_folder = get_report_s3_path(params)
    logging.info(f"report path: {report_folder}")

    # startup_script 
    startup_script = '''#!/bin/bash -e
trap 'catch $? $LINENO' EXIT
catch() {
  echo "catching!"
  if [ "$1" != "0" ]; then
    echo "Error $1 occurred on $2" > /tmp/ERROR
    '''
    startup_script += f"{params['aws_exe']} s3 cp /tmp/ERROR {report_folder}\n"
    startup_script += '''
  fi
}
git -C /tmp clone https://github.com/phenopolis/nextflow-genome;
'''

    what_to_do = set([i.strip() for i in params['do'].split(',')])
    if len(what_to_do) == 0:
        raise ValueError('Please specify what to do using --do')
    unknown_do = what_to_do - set(COMPONENTS.values())
    if len(unknown_do) > 0:
        raise ValueError(f'Some unknown components fed to --do "{unknown_do}". Only accept one or more of align, variantCall and annotation')
    if not params['input_table']:
        raise ValueError('Please specify --input_table (must be a S3 path)')

    to_do_list = [k for k,v in COMPONENTS.items() if v in what_to_do]
    for component_index in to_do_list:
        # make a s3 path to put sentinel file once pipeline is finished
        component = COMPONENTS[component_index]
        bucket_dir = f"{params['bucket_dir_base']}/{task_name}/{component}"
        startup_script = startup_script + f"cd /tmp/nextflow-genome/{component}; nextflow run {component}.nf --input_table {params['input_table']} -bucket-dir {bucket_dir} -with-report {report_folder}/{component}.html;"
        logging.info(f"bucket-dir for {component} will be: {bucket_dir}")

    startup_script += f"touch /tmp/DONE && {params['aws_exe']} s3 cp /tmp/DONE {report_folder}"

    # upload startup_script to s3
    s3_conn = boto3.client('s3')
    bare_report_path = report_folder.lstrip('s3://')
    bucket, prefix = bare_report_path.split('/', 1)
    local_tmp_path = os.path.join('tmp', prefix)
    logging.info(f"local tmp path: {local_tmp_path}")

    os.makedirs(local_tmp_path, exist_ok=True)
    startup_file = os.path.join(local_tmp_path, 'startup.sh')
    startup_s3_file = f"{prefix}/startup.sh"
    with open(startup_file, 'wt') as outf:
        outf.write(startup_script)
    s3_conn.upload_file(startup_file, bucket, startup_s3_file)

    # launch instance
    ec2 = boto3.resource('ec2')
    user_data = '''#!/bin/bash
'''
    user_data += f"su ec2-user -c '{params['aws_exe']} s3 cp s3://{bucket}/{startup_s3_file} /tmp/startup.sh' \n"
    #user_data += f"su ec2-user -c '{params['aws_exe']} s3 cp s3://phenopolis-nextflow/test/launch.py /tmp/startup.sh'"
    #user_data += 'touch /tmp/mememe'
    user_data += "su ec2-user -c 'bash /tmp/startup.sh'"
    instance = ec2.create_instances(
      ImageId=params['ec2']['image_id'],
      MinCount=1, MaxCount=1,
      KeyName=params['ec2']['key_name'],
      Placement={
        'AvailabilityZone': params['ec2']['availability_zone']
      },
      SecurityGroupIds=[params['ec2']['security_group']],
      UserData=user_data,
      InstanceType=params['ec2']['instance_type'],
      SubnetId=params['ec2']['subnet_id']
    )
    logging.info(f"job dispatched to {instance[0]}. You can log to it using KeyName {params['ec2']['key_name']}.")

    # keep an eye on report-folder
    wait = True
    with_error = False
    while wait:
        response = s3_conn.list_objects(Bucket=bucket, Prefix=prefix)
        if 'Content' in response:
            files = [i['Key'] for i in s3_conn.list_objects(Bucket=bucket, Prefix=prefix)['Contents']]
            files = set([i.split('/')[-1] for i in files])
            if files & set(['DONE', 'ERROR']):
                wait = False
                if 'ERROR' in files:
                    with_error = True
            else:
                time.sleep(params['check_interval'])
        else:
            time.sleep(params['check_interval'])

    logging.info(f"job done.")
    if with_error:
        logging.warning(f"job finishes with error!")
        

    logging.info(f"terminating instance")
    response = boto3.client.terminate_instances(
        InstanceIds=[
            instance[0].id,
        ],
    )
    logging.info(f"{response}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--do', dest='do',
        help='comma separated compoents of the pipeline you wish to run (e.g. --do align,variantCall,annotation')
    parser.add_argument('--input_table', dest='input_table',
        help='(S3) path to a tab-delimited file, withOUT header, with columns as (PID fastq1path fastq2path)')
    opts = parser.parse_args()
    with open('config.yml', 'rt') as inf:
        params = yaml.safe_load(inf)
    
    for key in ('do', 'input_table'):
        params[key] = getattr(opts, key)

    main(params)

