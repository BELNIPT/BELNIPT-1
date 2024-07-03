#!/usr/bin/env nextflow

// Script parameters
params.input = ""
params.output = ""

log.info """\
Input  : $params.input
Output : $params.output
"""

process wisecondor {

  conda '/conda/envs/wisecondor'
  publishDir "$params.output", mode: 'copy'

  input:
    path input from params.input

  output:
    path "reference_output.npz" into wisecondor_refset

  script:
  """
  #generate npz
  WisecondorX newref $input/*.npz reference_output.npz
  """

}

workflow.onComplete {
  """
  echo "DONE!"
  """
  println "Pipeline completed"
}
