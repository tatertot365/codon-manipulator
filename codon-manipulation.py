#!/usr/bin/env python 3
import sys
import re




def ReplaceSpace(sen, char):
    sentence = ""
    for i in sen:
        if i == char:
            continue
        else:
            sentence = sentence + i
    return sentence


def Length(var):
    j = 0
    for i in var:
        j += 1
    return j


def Alpha(sequence):
    seq1 = ""
    for i in sequence:
        if i == "a":
            seq1 = seq1 + "A"
        if i == "c":
            seq1 = seq1 + "C"
        if i == "t":
            seq1 = seq1 + "T"
        if i == "g":
            seq1 = seq1 + "G"
        if i not in ["a", "c", "g", "t"]:
            seq1 = seq1 + i
    return seq1


def calcGC (header, valid_sequence):
    key1 = ["G", "g"]
    key2 = ["C", "c"]
    GCount = 0
    CCount = 0
    for i in range(Length(valid_sequence)):
        if valid_sequence[i] in key1:
            GCount = GCount + 1
        if valid_sequence[i] in key2:
            CCount = CCount + 1
    GCCount = GCount + CCount
    percent = (GCCount / Length(valid_sequence)) * 100
    return header[1::] + "\t" + str(percent) + "%\n"


def revComp (header, valid_sequence):
    comp = ""
    for i in valid_sequence:
        if i == "A":
            comp = comp + "T"
        elif i == "C":
            comp = comp + "G"
        elif i == "G":
            comp = comp + "C"
        elif i == "T":
            comp = comp + "A"
    rev_comp = comp[::-1]
    rev_comp = "".join(rev_comp)
    return header + "\n" + rev_comp + "\n"


def transcribe(header, valid_sequence):
    valid_sequence = Alpha(valid_sequence)
    regex = r"T"
    replace = r"U"
    trans = re.sub(regex, replace, valid_sequence)
    return header + "\n" + trans + "\n"


def translate (header, valid_sequence, codons_dic):
    valid_sequence = Alpha(valid_sequence)
    regex = r"T"
    replace = r"U"
    trans = re.sub(regex, replace, valid_sequence)
    amino = ""
    protein = ""
    j = 0
    while j < Length(trans):
        amino = trans[j] + trans[j + 1] + trans[j + 2]
        if amino in codons_dic.keys():
            if codons_dic[amino] == "*":
                break
            else:
                protein = protein + codons_dic.get(amino)
        j += 3

    return header + "\n" + protein + "\n"


def getCount(str_to_look_in, char_to_count):
    char_count = 0
    for i in str_to_look_in:
        if i == char_to_count:
            char_count += 1
    return char_count


def countNucs (header, valid_sequence):
    length = Length(valid_sequence)
    a_count = getCount(valid_sequence, "A")
    c_count = getCount(valid_sequence, "C")
    t_count = getCount(valid_sequence, "T")
    g_count = getCount(valid_sequence, "G")
    perA = (a_count / length) * 100
    perC = (c_count / length) * 100
    perT = (t_count / length) * 100
    perG = (g_count / length) * 100
    return header[1::] + "\t" + str(length) + "\t" + str(a_count) + "(" + str(perA) + "%)\t" + str(c_count) + "(" + str(perC) + "%)\t" + str(
        g_count) + "(" + str(perG) + "%)\t" + str(t_count) + "(" + str(perT) + "%)\n"


def incorrect (numArgs):
    LS = " ".join(numArgs)
    print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
    print("\tAvailable options:")
    print("\t\t-g GC percent")
    print("\t\t-r reverse complement")
    print("\t\t-s transcription")
    print("\t\t-l translation")
    print("\t\t-n count nucleotides")
    sys.exit()


def errorCheckHead (line):
    write_line = True
    for i in line:
        replaced = ReplaceSpace(line, " ")
        for j in replaced[1:]:
            if j in ["A", "C", "G", "T"]:
                write_line = False
            else:
                write_line = True
                break
        if write_line:
            return True
        if not write_line:
            return False


Arguments = sys.argv[1:]

if Length(Arguments) == 5 or Length(Arguments) == 7:
    command = sys.argv[1]
    indash = sys.argv[2]
    infile = sys.argv[3]
    outdash = sys.argv[4]
    outfile = sys.argv[5]
    if indash == "-o":
        incorrect(Arguments)
    elif indash == "-i" and outdash == "-o" and command in ["-g", "-r", "-s", "-l", "-n", "-c"]:
        if command == "-g":
            with open(infile) as toRead, open(outfile, "w") as toWrite:
                toWrite.write("ID\tGC%\n")
                header = ""
                for line in toRead:
                    if line[0] == ">":
                        line = line.strip()
                        if errorCheckHead(line):
                            header = line
                        if not errorCheckHead(line):
                            toWrite.write(header[1::] + "\tERROR\n")
                    else:
                        error = False
                        upper = Alpha(line)
                        upper = upper.strip()
                        shotty = ""
                        for i in upper:
                            if i in ["A", "C", "T", "G"]:
                                shotty = shotty + i
                            elif i in [" ", "\t"]:
                                continue
                            else:
                                error = True
                        if error:
                            toWrite.write(header[1::] + "\tERROR" + "\n")
                        else:
                            write = calcGC(header, shotty)
                            toWrite.write(write)
        if command == "-r":
            with open(infile) as toRead, open(outfile, "w") as toWrite:
                header = ""
                for line in toRead:
                    if line[0] == ">":
                        line = line.strip()
                        if errorCheckHead(line):
                            header = line
                        if not errorCheckHead(line):
                            toWrite.write(header + "\nERROR\n")

                    else:
                        error = False
                        upper = Alpha(line)
                        upper = upper.strip()
                        shotty = ""
                        for i in upper:
                            if i in ["A", "C", "T", "G",]:
                                shotty = shotty + i
                            elif i in [" ", "\t"]:
                                continue
                            else:
                                error = True
                        if error:
                            toWrite.write(header + "\nERROR" + "\n")
                        else:
                            write = revComp(header, shotty)
                            toWrite.write(write)
        if command == "-s":
            with open(infile) as toRead, open(outfile, "w") as toWrite:
                header = ""
                for line in toRead:
                    if line[0] == ">":
                        line = line.strip()
                        if errorCheckHead(line):
                            header = line
                        if not errorCheckHead(line):
                            toWrite.write(header + "\nERROR\n")

                    else:
                        error = False
                        upper = Alpha(line)
                        upper = upper.strip()
                        shotty = ""
                        for i in upper:
                            if i in ["A", "C", "T", "G"]:
                                shotty = shotty + i
                            elif i in [" ", "\t"]:
                                continue
                            else:
                                error = True
                        if error:
                            toWrite.write(header + "\nERROR" + "\n")
                        else:
                            write = transcribe(header, shotty)
                            toWrite.write(write)
        if command == "-l":
            if Length(Arguments) == 7:
                amino_acid = {}
                with open(infile) as toRead, open(outfile, "w") as toWrite, open(sys.argv[7]) as codons:
                    for line in codons:
                        line = line.strip()
                        columns = line.split("\t")
                        codon = columns[0]
                        amino = columns[1]
                        amino_acid[codon] = amino
                    header = ""
                    for line in toRead:
                        if line[0] == ">":
                            line = line.strip()
                            if errorCheckHead(line):
                                header = line
                            if not errorCheckHead(line):
                                toWrite.write(header + "\nERROR\n")

                        else:
                            error = False
                            upper = Alpha(line)
                            upper = upper.strip()
                            shotty = ""
                            for i in upper:
                                if i in ["A", "C", "T", "G",]:
                                    shotty = shotty + i
                                elif i in [" ", "\t"]:
                                    continue
                                else:
                                    error = True
                            if error:
                                toWrite.write(header + "\nERROR" + "\n")
                            else:
                                write = translate(header, shotty, amino_acid)
                                toWrite.write(write)
            else:
                incorrect(Arguments)
        if command == "-n":
            with open(infile) as toRead, open(outfile, "w") as toWrite:
                toWrite.write("ID\tLength\tA(%A)\tC(%C)\tG(%G)\tT(%T)\n")
                header = ""
                for line in toRead:
                    if line[0] == ">":
                        line = line.strip()
                        if errorCheckHead(line):
                            header = line
                        if not errorCheckHead(line):
                            toWrite.write(header[1::] + "\tERROR\n")

                    else:
                        error = False
                        upper = Alpha(line)
                        upper = upper.strip()
                        shotty = ""
                        for i in upper:
                            if i in ["A", "C", "T", "G"]:
                                shotty = shotty + i
                            elif i in [" ", "\t"]:
                                continue
                            else:
                                error = True
                        if error:
                            toWrite.write(header[1::] + "\tERROR" + "\n")
                        else:
                            write = countNucs(header, shotty)
                            toWrite.write(write)
    else:
        incorrect(Arguments)
else:
    incorrect(Arguments)
