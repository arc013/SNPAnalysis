from flask import Flask, render_template, request, redirect, url_for, flash, session
import os
from werkzeug.utils import secure_filename
import sys
sys.path.append('scripts/')
from parseSnpedia import getrsIDfromVCFFile 
from snpedia import *
import ast
import subprocess
import math
import parseAnnotation
import myvariant_api

UPLOAD_FOLDER = 'uploads/'
ALLOWED_EXTENSIONS = set(['vcf'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config.update(dict(
    SECRET_KEY='sess'
))

snps = {}
summaries = {}

def get_files(path):
    filenames = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".csv"):
                filenames.append(os.path.join(root,file))   
    return filenames

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def main():
    if request.method == "POST":
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No file selected')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if request.form.getlist("check") and request.form.getlist('check')[0] == "on":
                toRun = ["java",'-Xmx4g',"-jar", "snpEff/snpEff.jar", "GRCh38.86", "uploads/" + file.filename, ">", "uploads/" + file.filename.rsplit('.', 1)[0] + ".annotation.vcf"]
                os.system(" ".join(toRun))
            return redirect(url_for('uploaded_file', filename=filename, index=1))
        else:
            flash("Invalid file extension")
            return redirect(request.url)
    else:
        return render_template('main.html')

def getSummary(rsid):
    res = myvariant_api.querey_data(rsid)
    toReturn = "No condition information available"
    if res == None:
        return toReturn
    if "disease_association" in res:
        assoc = res["disease_association"]
        if type(assoc) == list and assoc != []:
            toReturn = assoc[0]["condition"]
        elif "condition" in toReturn:
            toReturn = assoc["condition"]

    return toReturn

@app.route('/<filename>/<index>', methods=["GET", "POST"])
def uploaded_file(filename, index):
    if request.method == "POST":
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if request.form.getlist("check") and request.form.getlist('check')[0] == "on":
                toRun = ["java",'-Xmx4g',"-jar", "snpEff/snpEff.jar", "GRCh38.86", "uploads/" + file.filename, ">", "uploads/" + file.filename.rsplit('.', 1)[0] + ".annotation.vcf"]
                os.system(" ".join(toRun))
                if os.path.isfile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf"):
                    return redirect(url_for('uploaded_file_levelName', filename=filename, check_keyword="all_levels", index=1))
            return redirect(url_for('uploaded_file', filename=filename, index=1))
        else:
            flash("Invalid file extension")
            return redirect(request.url)
    if os.path.isfile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf"):
        return redirect(url_for('uploaded_file_levelName', filename=filename, check_keyword="all_levels", index=1))
    chrNums = getrsIDfromVCFFile("uploads/"+filename)
    rsidsDict = {}
    for chromosome, rsids in chrNums.items():
        for rsid in rsids:
            rsidsDict[(rsid, chromosome)] = True
    rsids = rsidsDict.keys()
    
    numPages = int(math.ceil(len(rsids)/50))
    currRsids = rsids[int(index)*50:int(index)*50+50]
    currRsids = [[rsid[0], rsid[1], getSummary(rsid[0])] for rsid in currRsids]
    canNext = int(index) < numPages
    canNext2 = int(index) < numPages-1
    canPrev = int(index) > 0
    canPrev2 = int(index) > 1
    canFilter = False
    if os.path.isfile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf"):
        canFilter = True
    return render_template("view.html", canFilter=canFilter, index=int(index), filename=filename, rsids=currRsids, numPages=numPages, canNext=canNext, canPrev=canPrev, canPrev2=canPrev2, canNext2=canNext2)

@app.route('/name/<filename>/<check_keyword>/<index>', methods=["GET", "POST"])
def uploaded_file_levelName(filename, check_keyword, index):
    if request.method == "POST":
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if request.form.getlist("check") and request.form.getlist('check')[0] == "on":
                toRun = ["java",'-Xmx4g',"-jar", "snpEff/snpEff.jar", "GRCh38.86", "uploads/" + file.filename, ">", "uploads/" + file.filename.rsplit('.', 1)[0] + ".annotation.vcf"]
                os.system(" ".join(toRun))
                if os.path.isfile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf"):
                    return redirect(url_for('uploaded_file_levelName', filename=filename, check_keyword="all_levels", index=1))
            return redirect(url_for('uploaded_file', filename=filename, index=1))
        else:
            flash("Invalid file extension")
            return redirect(request.url)
    chrNums = getrsIDfromVCFFile("uploads/"+filename)
    every=False
    origCheckKeyword = check_keyword
    if(check_keyword == "all_levels"):
        check_keyword = ""
        every = True
    rsids = parseAnnotation.getrsIDfromVCFFile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf", "order_annotations.txt", check_keyword=check_keyword, get_pop=True, every=every)
    if(rsids.keys() == []):
        flash(filename + " does not contain the keyword " + check_keyword)
        return redirect(url_for('uploaded_file', filename=filename, index=index))
    rsResList = []
    for chromosome, val in rsids.items():
        for rsTuple in val:
            rsRes = []
            rsRes.append(rsTuple[0])
            rsRes.append(None)
            rsRes.append(rsTuple[1])
            rsRes.append(chromosome)
            if "EUR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["EUR_AF"])
            else:
                rsRes.append(0)
            if "SAS_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["SAS_AF"])
            else:
                rsRes.append(0)
            if "AFR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["AFR_AF"])
            else:
                rsRes.append(0)
            if "AMR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["AMR_AF"])
            else:
                rsRes.append(0)
            if "EAS_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["EAS_AF"])
            else:
                rsRes.append(0)
            rsResList.append(rsRes)
    rsids = rsResList
    numPages = int(math.ceil(len(rsids)/50))
    currRsids = rsids[int(index)*50:int(index)*50+50]
    for i, res in enumerate(currRsids):
        currRsids[i][1] = getSummary(res[0].title())
    canNext = int(index) < numPages
    canNext2 = int(index) < numPages-1
    canPrev = int(index) > 0
    canPrev2 = int(index) > 1
    return render_template("view_level.html", check_keyword=origCheckKeyword, type="_levelName", index=int(index), filename=filename, rsids=currRsids, numPages=numPages, canNext=canNext, canPrev=canPrev, canPrev2=canPrev2, canNext2=canNext2)

@app.route('/num/<filename>/<order_level>/<index>', methods=["GET", "POST"])
def uploaded_file_levelNum(filename, order_level, index):
    if request.method == "POST":
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if request.form.getlist("check") and request.form.getlist('check')[0] == "on":
                toRun = ["java",'-Xmx4g',"-jar", "snpEff/snpEff.jar", "GRCh38.86", "uploads/" + file.filename, ">", "uploads/" + file.filename.rsplit('.', 1)[0] + ".annotation.vcf"]
                os.system(" ".join(toRun))
                if os.path.isfile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf"):
                    return redirect(url_for('uploaded_file_levelName', filename=filename, check_keyword="all_levels", index=1))
            return redirect(url_for('uploaded_file', filename=filename, index=1))
        else:
            flash("Invalid file extension")
            return redirect(request.url)
    chrNums = getrsIDfromVCFFile("uploads/"+filename)
    if(order_level == "none"):
        order_level = 46
    rsids = parseAnnotation.getrsIDfromVCFFile("uploads/" + filename.rsplit('.', 1)[0] + ".annotation.vcf", "order_annotations.txt", order_level=int(order_level), get_pop=True)
    if(rsids.keys() == []):
        flash(filename + " does not contain the order level " + order_level)
        return redirect(url_for('uploaded_file', filename=filename, index=index))
    rsResList = []
    for chromosome, val in rsids.items():
        for rsTuple in val:
            rsRes = []
            rsRes.append(rsTuple[0])
            rsRes.append(None)
            rsRes.append(rsTuple[1])
            rsRes.append(chromosome)
            if "EUR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["EUR_AF"])
            else:
                rsRes.append(0)
            if "SAS_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["SAS_AF"])
            else:
                rsRes.append(0)
            if "AFR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["AFR_AF"])
            else:
                rsRes.append(0)
            if "AMR_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["AMR_AF"])
            else:
                rsRes.append(0)
            if "EAS_AF" in rsTuple[2]:
                rsRes.append(rsTuple[2]["EAS_AF"])
            else:
                rsRes.append(0)
            rsResList.append(rsRes)
    rsids = rsResList
    
    numPages = int(math.ceil(len(rsids)/50))
    currRsids = rsids[int(index)*50:int(index)*50+50]
    for i, res in enumerate(currRsids):
        currRsids[i][1] = getSummary(res[0].title())

    canNext = int(index) < numPages
    canNext2 = int(index) < numPages-1
    canPrev = int(index) > 0
    canPrev2 = int(index) > 1
    return render_template("view_level.html", order_level=order_level, type="_levelNum", index=int(index), filename=filename, rsids=currRsids, numPages=numPages, canNext=canNext, canPrev=canPrev, canPrev2=canPrev2, canNext2=canNext2)


@app.route('/<rsid>', methods=["GET", "POST"])
def getPubs(rsid):
    if request.method == "POST":
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if request.form.getlist("check") and request.form.getlist('check')[0] == "on":
                toRun = ["java",'-Xmx4g',"-jar", "snpEff/snpEff.jar", "GRCh38.86", "uploads/" + file.filename, ">", "uploads/" + file.filename.rsplit('.', 1)[0] + ".annotation.vcf"]
                os.system(" ".join(toRun))
            return redirect(url_for('uploaded_file', filename=filename, index=1))
        else:
            flash("Invalid file extension")
            return redirect(request.url)
    pubs = getSnpedia(rsid)["pubs"]
    dictionary = myvariant_api.querey_data("rs17822931")
    if "evs_clinical_info" in dictionary:
        pubs = dictionary["evs_clinical_info"].split("|")
        del dictionary["evs_clinical_info"]
    assoc = []
    if "disease_association" in dictionary:
        assoc = dictionary["disease_association"]
        for i, disease in enumerate(assoc):
            if "identifiers" in disease:
                assoc[i]["identifiers"] = ", ".join([key+"="+val for key, val in disease["identifiers"].items()])
        del dictionary["disease_association"]
    for key, val in dictionary.items():
        if type(val) == dict:
            dictionary[key] = ", ".join([key2+"="+val2 for key2, val2 in val.items()])
    summary = getSummary(rsid)
    snpRes = getSnpedia(rsid)
    if rsid.title() in summaries:
        if summaries[rsid.title()] != []:
            summary = summaries[rsid.title()][0]
    if snpRes and "genos" in snpRes and snpRes["genos"] != [] and "summary" in snpRes["genos"][0]:
        summary = snpRes["genos"][0]["summary"]
    pmed = []
    if "publication_pmid" in dictionary:
        pmed = dictionary["publication_pmid"]
        del dictionary["publication_pmid"]
    return render_template("pubs.html", rsid=rsid, pubs=pubs, d=dictionary, assoc=assoc, summary=summary, pmed=pmed)

def readSnps():
    list_files = get_files('all_rsids/')
    for filename in list_files:
        f = open(filename, "rb")
        for rsid in f.readlines():
            snps[rsid.strip()] = True
        f.close()


def readSummaries():
    f = open("summaries.csv", "rb")
    summariesCsv = csv.DictReader(f)
    for row in summariesCsv:
        summaries[row["genotype"]] = ast.literal_eval(row["summary"])
    f.close()
    
if __name__ == "__main__":
    readSnps()
    readSummaries()
    app.run(threaded=True)