import sys
import os
import logging
import LicenseUtils

logging.basicConfig(format="[UpdateLicense]: %(message)s", level = logging.INFO)
for filename in sys.stdin.readlines():
    filename = filename.strip()

    if LicenseUtils.isSourceFile(filename) and os.path.exists(filename):
        licenseFile = open(LicenseUtils.getAppropriateLicense(filename))
        sourceFile = open(filename)
        updatedSource = "/*\n"
        for line in licenseFile.readlines():
            updatedSource += "* " + line

        skipLicense = True
        for line in sourceFile.readlines():
            strippedLine = line.strip()
            if skipLicense and strippedLine.startswith("package"):
                updatedSource += "*/\n\n"
                skipLicense = False
            elif skipLicense and LicenseUtils.lineIsNotCommentedOut(strippedLine):
                logging.error(filename + " is missing package information")
                exit(2)

            if not skipLicense:
                updatedSource += line

        sourceFile.close()
        sourceFile = open(filename, "w")
        sourceFile.write(updatedSource)
        sourceFile.close()

        logging.info(filename + " [license successfully added]")
    else:
        logging.debug(filename + " doesn't require a license")
