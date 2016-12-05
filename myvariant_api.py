import sys
import requests
#print ("hi")

def querey_data(rs_id):
    print ("input rs_id is : "+rs_id)
    api_id = "http://myvariant.info/v1/query?q="+rs_id
    result = requests.get(api_id)
    hits = result.json()["hits"] #an array
    index = -1
    grasp_i = -1
    cadd_i  = -1
    evs_i   = -1
    for i in range (len(hits)):
        if "clinvar" in hits[i]:
            index = i
        if "grasp" in hits[i]:
            grasp_i = i
        if "cadd" in hits[i]:
            cadd_i = i
        if "evs" in hits[i]:
            evs_i = i
    if index == -1:
        return 
 
    jake = {}
    clinvar_object = hits[index]["clinvar"]
    if "hg38" in clinvar_object:
        start_end = str(clinvar_object["hg38"]["start"]) +"-"+str(clinvar_object["hg38"]["end"])
        jake["start-end"] = start_end
    if "variant_id" in clinvar_object:
        jake["variant_id"] = clinvar_object["variant_id"]
    if "cytogenic" in clinvar_object:
        jake["cytogenic"] = clinvar_object["cytogenic"]
    if "allele_id" in clinvar_object:
        jake["allele_id"] = clinvar_object["allele_id"]
    if "rcv" in clinvar_object:
        rcv_array = clinvar_object["rcv"]
        for i in rcv_array:
            jake["disease_association"]=[]
            condition = {}
            if "conditions" in i:
                if "name" in i["conditions"]:
                    condition["condition"]=i["conditions"]["name"]
                if "identifiers" in i["conditions"]:
                    condition["identifiers"] =i["conditions"]["identifiers"]
            if "clinical_significance" in i:
                condition["clinical_significance"]=i["clinical_significance"]
            if "accession" in i:
                condition["accession"] = i["accession"]
            jake["disease_association"].append(condition)
    if grasp_i != -1:
        grasp_object = hits[grasp_i]["grasp"]["publication"]
        pmid_array = []
        for i in grasp_object:
            pmid_array.append(i["pmid"])
        jake["publication_pmid"]=pmid_array
    if cadd_i != -1:
        if "gene" in hits[cadd_i]:
            cadd_gene = {}
            for i in hits[cadd_i]["gene"]:
                if "feature_id" in i:
                    cadd_gene["feature_id"]=i["feature_id"]
                if "gene_id" in i:
                    cadd_gene["gene_id"]=i["gene_id"]
                if "genename" in i:
                    cadd_gene["genename"]=i["genename"]
            jake["cadd_gene"] = cadd_gene

    if evs_i != -1:
        if "clinical_info" in hits[evs_i]["evs"]:
            jake["evs_clinical_info"] =  hits[evs_i]["evs"]["clinical_info"]
                
    print (jake)
    return jake


#querey_data()
    
 # hd38 {start-end: "123-456"}
# end and start benign, 
# variant_id,
#cyogenic,
#allele_id,

# clinical_significance,, accession, name, identifiers, 