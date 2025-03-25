import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { MSNBASXIC } from '/home/proteomics/mygit/nf-core/nf-core-ribomsqc/modules/local/msnbasexic/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
      input[0] = [
        [ id: 'test' ],
        file("/home/proteomics/mygit/nf-core/nf-core-ribomsqc/modules/local/msnbasexic/tests/data/input/20240315_QCN1_001_03_New_STD.mzML",checkIfExists: true),
        file("/home/proteomics/mygit/nf-core/nf-core-ribomsqc/modules/local/msnbasexic/tests/data/input/qcn1.tsv",checkIfExists: true)
      ]
      input[1] =  "m3C"
      input[2] =  150
      input[3] =  20
      input[4] =  2
      input[5] =  false
      input[6] =  false
      input[7] =  "xic_plot.pdf"
      input[8] =  true
      
    //----

    //run process
    MSNBASXIC(*input)

    if (MSNBASXIC.output){

        // consumes all named output channels and stores items in a json file
        for (def name in MSNBASXIC.out.getNames()) {
            serializeChannel(name, MSNBASXIC.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = MSNBASXIC.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
