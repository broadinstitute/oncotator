import re
import logging

def lineIsNotCommentedOut(line):
    return not line == "" and not line.startswith(" *") and not line.startswith("*") and not line.startswith("/*") and not line.startswith("//") and not line == "\n" and not line == "\r\n"

def isSourceFile(filename):
    return str.endswith(filename, ".java") or str.endswith(filename, ".scala")

def getAppropriateLicense(filename):
    licenseFileName = "broad_license_oncotator_oday.txt"
    return licenseFileName

def blankLine(line):
    return re.match("^$", line) is not None

def lineIsCommentedOut(line):
    return line.startswith("//") or line.startswith("/*") or line.startswith("*")

def extractLicenseFromSource(file):
    extractedLicense = ""
    inCommentBlock = 0
    for line in file.readlines():
        strippedLine = line.strip()
        if strippedLine.startswith("package"):
            return extractedLicense # found package line, return the license.

        if blankLine(strippedLine):
            continue                # just skip blank lines

        if strippedLine.startswith('"""'):
            inCommentBlock += 1     # mark start of a comment block

        if inCommentBlock or lineIsCommentedOut(strippedLine):
            extractedLicense += strippedLine.strip(" * ").strip(" // ")
        else:
            break

        if strippedLine.endswith("*/"):
            inCommentBlock -= 1     # mark end of a comment block

    return False

def extractLicense(filename):
    licenseFile = open(filename)
    license = ""
    for line in licenseFile.readlines():
        license += line.strip()
    return license

