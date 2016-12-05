from bs4 import BeautifulSoup
import urllib
import json
import sys
import csv
import codecs
import time

startLink = "http://bots.snpedia.com/index.php?title=Category:In_dbSNP&pagefrom="
endLink = "#mw-pages"

def atEnd(soup):
	seenNext = False
	for mwPage in soup.find_all("div", id="mw-pages"):
		for childA in mwPage.find_all('a'):
			if childA.text == "next page":
				seenNext = True
	return seenNext

def getRsIds(prevRs, rsIds):
	r = urllib.urlopen(startLink + prevRs + endLink).read()
	soup = BeautifulSoup(r, "html.parser")
	mwCategory = soup.find_all("div", class_="mw-category")
	rows = (mwCategory[0]).find_all("li")
	for row in rows:
		if row.text[:2] != "Rs":
			continue
		rsIds.append(row.text)
	canContinue = atEnd(soup)
	return canContinue, rsIds

def writeToFile(rsIds):
	f = open("all_rsids/"+rsIds[0] + ".csv", "wb")
	for rsId in rsIds:
		f.write("%s\n" % rsId)
	f.close()

prevRs = "Rs10"
rsIds = []
rsCount = 0

while True:
	rsCount += 1
	canContinue, rsIds = getRsIds(prevRs, rsIds)
	prevRs = "Rs" + str(int(rsIds[-1][2:])+1)
	if(rsCount % 500 == 499):
		writeToFile(rsIds)
		rsIds = []
	time.sleep(1)
	if canContinue == False:
		break
writeToFile(rsIds)