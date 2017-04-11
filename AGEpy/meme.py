def filterMotifs(memeFile,outFile, minSites):
    """
    Selectes motifs from a meme file based on the number of sites.

    :param memeFile: MEME file to be read
    :param outFile: MEME file to be written
    :param minSites: minimum number of sites each motif needs to have to be valid

    :returns: nothing
    """

    with open(memeFile, "r") as mF:
        oldMEME=mF.readlines()
        newMEME=oldMEME[:7]
        i=7
        while i < len(oldMEME):
            if oldMEME[i].split(" ")[0] == "MOTIF":
                print oldMEME[i].split("\n")[0], int(oldMEME[i+2].split("nsites= ")[1].split(" ")[0])
                sys.stdout.flush()
                if int(oldMEME[i+2].split("nsites= ")[1].split(" ")[0]) > minSites:
                    newMEME.append(oldMEME[i])
                    f=i+1
                    while oldMEME[f].split(" ")[0] != "MOTIF":
                        newMEME.append(oldMEME[f])
                        f=f+1
                    i=i+1
                else:
                    i=i+1
            else:
                i=i+1
    with open(outFile, "w+") as out:
        out.write("".join(newMEME) )

    return newMEME
