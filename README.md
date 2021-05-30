# Use Nextflow to dispatch Exome/Genome pipelines to AWS Batch

Add aws credentials before running the script as:
`nextflow run main.nf [-resume] [-with-report report.html]`

![deployment_options](https://user-images.githubusercontent.com/47983965/119678496-7926d980-be37-11eb-95b5-a1d8db62e0e0.png)

![aws_deployment](https://user-images.githubusercontent.com/47983965/119678501-7a580680-be37-11eb-9518-0919cbf72340.png)

# Arrangement of the pipeline
The Pipeline is split into three parts: `align`, `variantCall` and `annotation`. User can also choose to run `exomiser` following `annotation`

To run, go to individal folder and run `nextflow run [pipeline name].nf`

# Setup of an API server
The flask folder has a script to lunch an API server to accept tasks either to run the complete set of the pipeline, or part of it.

To run, go to the flask folder and run `python main.py` (probably want to edit the `config.yml` in the same folder to adapt your working evironment).
