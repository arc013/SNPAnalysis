from bs4 import BeautifulSoup
import urllib
import json
import sys
import csv
import codecs

def getSummaryInfo(offset, summaries):
	url = "http://www.snpedia.com/index.php/Special:Ask/-5B-5BCategory:Is-20a-20genotype-5D-5D/-3FMagnitude/-3FRepute/-3FSummary/mainlabel%3D/offset%3D" + str(offset) + "/limit%3D500/order%3DDESC/sort%3D/prettyprint%3Dtrue/format%3Djson"
	r = urllib.urlopen(url).read()
	r = json.loads(r)
	for genotype, info in r["results"].items():
		if genotype[0] != "R":
			continue
		summaries.append({"genotype": genotype.split("(")[0], "summary": info["printouts"]["Summary"]})
	return summaries
summaries = []
for offset in range(0, 5500, 500):
	print("On: " + str(offset))
	summaries = getSummaryInfo(offset, summaries)

print(len(summaries))
f = open("summaries.csv", "wb")
csvWriter = csv.DictWriter(f, fieldnames=["genotype", "summary"])
csvWriter.writeheader()
csvWriter.writerows(summaries)
f.close()
