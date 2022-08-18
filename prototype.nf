nextflow.enable.dsl=2


process pythonTest
{
    output:
    stdout

    script:
    """
    $scripts_folder/test.py
    """
}

process rTest
{
    output:
    stdout

    script:
    """
    $scripts_folder/test.py
    """
}

workflow  
{
    pythonTest().view()
    rTest().view()
}