import os
import sys
from datetime import datetime
import subprocess
from flask import Flask, request
from flask_cors import CORS, cross_origin
import yaml
import uuid
from threading import Thread

with open('config.yml', 'rt') as inf:
    params = yaml.safe_load(inf)

app = Flask(__name__)
CORS(app)

COMPONENTS = {
    1: 'align',
    2: 'variantCall',
    3: 'annotation',
}


weblogs = {}
def launcher(params):
    cwd = os.path.dirname(os.path.realpath(__file__))
    what_to_do = set([i.strip() for i in params['do'].split(',')])
    now = datetime.today().strftime('%Y-%m-%d-%H-%M-%S')
    report_dir = os.path.join('report', now)
    os.makedirs(report_dir, exist_ok=True)
    for do_what in ('align', 'variantCall', 'annotation'):
        if do_what in what_to_do:
            wd = os.path.join(cwd, '..', do_what)
            cmd = ['nextflow', '-C', os.path.join(wd, 'nextflow.config'), '-log', f"{report_dir}/nextflow-{do_what}.log", 'run', os.path.join(wd, f'{do_what}.nf'), '--input_table', params['input_table'], '-bucket-dir', f"{params['bucket_dir_base']}/{now}/{do_what}", '-with-report', f"{report_dir}/ report-{do_what}.html", '-with-weblog', f"{params['web_log_url']}/{params['job_id']}"]
            subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr).communicate()

@app.route("/run", methods = ['POST'])
@cross_origin()
def launch():
    data = request.form
    if not {'do', 'input_table'}.issubset(set(data.keys())):
        # needs better error handle
        return 'Error: needs both do: [align, variantCall, annotation] and input_table'
    params.update(data)
    params['job_id'] = uuid.uuid4()
    Thread(target=launcher, args=(params, )).start()
    return {'job_id':params['job_id']}

@app.route("/weblog/<job_id>", methods = ['GET','POST'])
@cross_origin()
def weblog(job_id):
    if request.method == 'GET':
        if job_id not in weblogs:
            return {'status': 'fail', 'message': 'No just job_id'}
        else:
            return {'status': 'success', 'data':weblogs[job_id]}
    if request.method == 'POST':
        weblogs[job_id] = request.json
        print('form')
        print(request.form)
        print('json')
        print(request.json)
        print('args')
        print(request.args)

        return {'status': 'success', 'message': 'data received'}

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=params['port'])
