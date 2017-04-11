
def getFasta(opened_file, sequence_name):
    """
    Retrieves a sequence from an opened multifasta file

    :param opened_file: an opened multifasta file eg. opened_file=open("/path/to/file.fa",'r+')
    :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)

    returns: a string with the sequence of interest
    """

    lines = opened_file.readlines()
    seq=str("")
    for i in range(0, len(lines)):
        line = lines[i]
        if line[0] == ">":
            fChr=line.split(" ")[0].split("\n")[0]
            fChr=fChr[1:]
            if fChr == sequence_name:
                s=i
                code=['N','A','C','T','G']
                firstbase=lines[s+1][0]
                while firstbase in code:
                    s=s + 1
                    seq=seq+lines[s]
                    firstbase=lines[s+1][0]

    if len(seq)==0:
        seq=None
    else:
        seq=seq.split("\n")
        seq="".join(seq)

    return seq

def writeFasta(sequence, sequence_name, output_file):
    """
    Writes a fasta sequence into a file.

    :param sequence: a string with the sequence to be written
    :param sequence_name: name of the the fasta sequence
    :param output_file: /path/to/file.fa to be written

    :returns: nothing
    """
    i=0
    f=open(output_file,'w')
    f.write(">"+str(sequence_name)+"\n")
    while i <= len(sequence):
        f.write(sequence[i:i+60]+"\n")
        i=i+60
    f.close()

def rewriteFasta(sequence, sequence_name, fasta_in, fasta_out):
    """
    Rewrites a specific sequence in a multifasta file while keeping the sequence header.

    :param sequence: a string with the sequence to be written
    :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
    :param fasta_in: /path/to/original.fa
    :param fasta_out: /path/to/destination.fa

    :returns: nothing
    """
    f=open(fasta_in, 'r+')
    f2=open(fasta_out,'w')
    lines = f.readlines()
    i=0
    while i < len(lines):
        line = lines[i]
        if line[0] == ">":
            f2.write(line)
            fChr=line.split(" ")[0]
            fChr=fChr[1:]
            if fChr == sequence_name:
                code=['N','A','C','T','G']
                firstbase=lines[i+1][0]
                while firstbase in code:
                    i=i+1
                    firstbase=lines[i][0]
                s=0
                while s <= len(sequence):
                    f2.write(sequence[s:s+60]+"\n")
                    s=s+60
            else:
                i=i+1
        else:
            f2.write(line)
            i=i+1

    f2.close
    f.close
