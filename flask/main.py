import os
import sys
from datetime import datetime
import subprocess
from flask import Flask, request
from flask_cors import CORS, cross_origin
import requests
import yaml
import uuid
from threading import Thread

with open('config.yml', 'rt') as inf:
    params = yaml.safe_load(inf)

app = Flask(__name__)
CORS(app)

COMPONENTS = (
    'align',
    'variantCall',
    'annotation',
    'exomiser',
)


weblogs = {}
def launcher(params):
    cwd = os.path.dirname(os.path.realpath(__file__))
    what_to_do = set([i.strip() for i in params['do'].split(',')])
    now = datetime.today().strftime('%Y-%m-%d-%H-%M-%S')
    report_dir = os.path.join(params['report_dir_base'], now)
    os.makedirs(report_dir, exist_ok=True)
    for do_what in COMPONENTS:
        if do_what in what_to_do:
            print(f"do {do_what}")
            wd = os.path.join(cwd, '..', do_what)
            cmd = ['nextflow', '-C', os.path.join(wd, 'nextflow.config'), '-log', f"{report_dir}/nextflow-{do_what}.log", 'run', os.path.join(wd, f'{do_what}.nf'), '--input_table', params['input_table'], '-bucket-dir', f"{params['bucket_dir_base']}/{do_what}", '-with-report', f"{report_dir}/ report-{do_what}.html", '-with-weblog', f"{params['web_log_url']}/{params['job_id']}", "-resume"]
            subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr).communicate()
            print(f"{do_what} complete")
    requests.post(f"{params['web_log_url']}/{params['job_id']}", json={'status': 'done', 'message': f"reports can be found at {report_dir}"})

@app.route("/run", methods = ['POST'])
@cross_origin()
def launch():
    data = request.form
    if not {'do', 'input_table'}.issubset(set(data.keys())):
        # needs better error handle
        return f"Error: needs both do: [{', '.join(COMPONENTS)}] and input_table"
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
        elif weblogs[job_id].get('status', None) == 'done':
            return {'status': 'done', 'data': weblogs[job_id]}
        else:
            return {'status': 'running', 'data':weblogs[job_id]}
    if request.method == 'POST':
        weblogs[job_id] = request.json

        return {'status': 'success', 'message': 'data received'}

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=params['port'])
