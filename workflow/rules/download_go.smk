rule download_go:
    output:
        go_json=config["go_json"],
        go_obo=config["go_obo"],
        gaf=config["gaf"]
    shell:
        """
        wget -qO {output.go_json} "http://purl.obolibrary.org/obo/go/go-basic.json"
        wget -qO {output.go_obo} "http://purl.obolibrary.org/obo/go/go-basic.obo"
        wget -qO {output.gaf}.gz "https://current.geneontology.org/annotations/goa_human.gaf.gz"
        gunzip -c {output.gaf}.gz > {output.gaf} && rm {output.gaf}.gz
        """
