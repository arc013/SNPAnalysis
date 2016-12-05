from bs4 import BeautifulSoup
import urllib
import json
import sys
import csv
import codecs

def writeFile(fileName, fieldNames, data):
	f = open(fileName, "wb")
	csvWriter = csv.DictWriter(f, fieldnames=fieldNames)
	csvWriter.writeheader()
	csvWriter.writerows(data)
	f.close()

def getSnpedia(geneId):
	r = urllib.urlopen('http://www.snpedia.com/index.php/'+geneId).read()
	soup = BeautifulSoup(r, "html.parser")
	return {"genos": getGenos(soup), "pubs": getPubs(soup)}

def getGenos(soup):
	smwtable = soup.find_all("table", class_="smwtable")
	if smwtable == []:
		return []
	rows = (smwtable[0]).find_all("tr")
	genos = []
	for row in rows:
		cells = row.find_all("td")
		if(cells == []):
			continue
		geno = []
		for cell in cells:
			geno.append(cell.text.strip())
		genos.append(geno)
	genos = [{"geno": geno[0], "mag": geno[1], "summary": geno[2]} for geno in genos]
	return genos

def getPubs(soup):
	content = soup.find_all("div", id="mw-content-text")
	ps = content[0].find_all("p")
	pubIds = []
	for p in ps:
		hrefs = p.find_all("a", class_="mw-magiclink-pmid")
		if hrefs == []:
			continue
		href = hrefs[0]["href"][2:]
		pubId = href.split("?")[0].split("/")[-1]
		pubIds.append(pubId)
	pubsStr = ",".join(pubIds)
	pubs = urllib.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + pubsStr).read()
	pubs = json.loads(pubs)
	if "result" not in pubs:
		return []
	pubsResult = pubs["result"]
	pubs = [{"id": pubId, "title": pubsResult[pubId]["title"].encode('ascii','ignore')} for pubId in pubIds]
	return pubs

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("You forgot to pass in a gene id...")
	geneId = sys.argv[1]
	r = urllib.urlopen('http://www.snpedia.com/index.php/'+geneId).read()
	soup = BeautifulSoup(r, "html.parser")
	genos = getGenos(soup)
	writeFile("genes/"+geneId+"_genos.csv", ["geno", "mag", "summary"], genos)
	pubs = getPubs(soup, True)
	writeFile("genes/"+geneId+"_pubs.csv", "wb", "utf-8", ["id", "title"], pubs)
