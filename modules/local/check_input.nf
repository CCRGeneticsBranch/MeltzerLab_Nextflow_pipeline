process Check_Input {

input:
path(samplesheet)

output:
path("samplesheet_checked.csv")

script:

"""

check_samplesheet.py ${samplesheet} samplesheet_checked.csv

"""

}
