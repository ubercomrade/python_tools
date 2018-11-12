import math
import csv
import sys
import os
import random
import itertools
import argparse

def readFasta(path, everyStr = False):
    '''
    Чтение мотивов из фаила, если это просто список строк с разделителем \n, то everyStr = True,
    если это fasta формат, то everyStr = False
    '''
    sequences = []
    if everyStr:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file]
    else:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file if i.strip() != '>']
    return(sequences)


def removeEqualentSeq(seqList, homology=0.95):
    '''
    Удаление гомологичных последовательностей из списка (seqList)
    Если кол-во совпадений при сравнении последовательности 1 и 2 >= длина последовательности * homology, 
    то последовательность 1 удаляется из списка
    Функция возвращает новый список
    '''
    seqList = list(seqList)
    treshold = homology * len(seqList[0])
    for seq1 in tuple(seqList):
        sub_seqList = list(seqList)
        sub_seqList.remove(seq1)
        for seq2 in sub_seqList:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seqList.remove(seq1)
                break
    return(seqList)


def removeEqualentSeq_2(seqList, homology=0.95):
    '''
    Медленный аналог функции removeEqualentSeq
    '''
    allMotifs = list(seqList)
    treshold = homology * len(allMotifs[0])
    for i in tuple(allMotifs):
        subMotifs = list(allMotifs)
        subMotifs.remove(i)
        for j in subMotifs:
            score = int()
            for l in range(len(i)):
                if score >= treshold:
                    break
                if i[l] == j[l]:
                    score += 1
            if score >= treshold:
                allMotifs.remove(i)
                break
    return(allMotifs)


def makePFMfromPCM(PCM, kind, pseudocount= '1/N'):
    '''
    Вычисление частотной матрицы на основе PCM.
    Для того чтобы избавиться от 0 значений частот используется pseudocount.
    Pseudocount может быть dict со стандартными значениями {'A':0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} [1],
    либо pseudocount может быть str со значением sqroot [2].
    Подробнее о расчетах смотри:
    1)Wyeth W.Wasserman and Albin Sandelin 
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315
      
    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a 
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481
      
    В любых других условиях функция ничего не возвращает
    
    '''
    
    NumberOfSites = int()
    for i in PCM.keys():
        NumberOfSites += PCM[i][0]
    
    if kind == 'di':
        PFM = {}
        diNucleotides = itertools.product('ACGT', repeat=2)
        for i in diNucleotides:
            PFM[''.join(i)] = []
    elif kind == 'mono':
        PFM = {}
        monoNucleotides = itertools.product('ACGT', repeat=1)
        for i in monoNucleotides:
            PFM[i[0]] = []
    else:
        print('ALARM!')
        pass
    
    if pseudocount == '1/N':
        firstKey = list(PCM.keys())[0]
        nucPseudo = 1/len(PCM.keys())
        for i in range(len(PCM[firstKey])):
            for nuc in PCM.keys():
                PFM[nuc].append((PCM[nuc][i] + nucPseudo) / (NumberOfSites + 1))            
        return(PFM)
    
    elif pseudocount == 'sqroot':
        totalSqRoot = int()
        for i in PCM.keys():
            totalSqRoot += PCM[i][0]
        totalSqRoot = math.sqrt(totalSqRoot)
        sqRoot = totalSqRoot/len(PCM.keys())
        
        firstKey = list(PCM.keys())[0]
        for i in range(len(PCM[firstKey])):
            for nuc in PCM.keys():
                PFM[nuc].append((PCM[nuc][i] + sqRoot) / (NumberOfSites + totalSqRoot))
            
        return(PFM)
    else:
        print('ALARM!')
        pass


def makePWMfromPCM(PCM, kind, background, method='log-odds', pseudocount='1/N'):
    '''
    Функиця, которая считает PWM (position weight matrix) на основе PCM (position count matrix)
    с преобразованием log-odds (добавить новые)
    
    Ref:
    1)Wyeth W.Wasserman and Albin Sandelin 
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315
      
    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a 
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481
      
    3)Oliver D. King and Frederick P. Roth
      A non-parametric model for transcription factor binding sites
      doi: 10.1093/nar/gng117
    
    '''
    if kind == 'di':
        PWM = {}
        diNucleotides = itertools.product('ACGT', repeat=2)
        for i in diNucleotides:
            PWM[''.join(i)] = []
    elif kind == 'mono':
        PWM = {}
        monoNucleotides = itertools.product('ACGT', repeat=1)
        for i in monoNucleotides:
            PWM[i[0]] = []
    else:
        print('ALARM!')
        pass
    
    PFM = makePFMfromPCM(PCM, kind, pseudocount)
    firstKey = list(PCM.keys())[0]
    for i in range(len(PFM[firstKey])):
        for j in PFM.keys():
            PWM[j].append(math.log(PFM[j][i] / background[j]))        
    return(PWM)


def makePCM(motifs, kind):
    '''
    input - список мотивов одинаковой длины
    output -  PCM
    kind is type of matrix di or mono
    
    Создает PCM на основе списка мотивов
    '''
    if kind == 'di':
        matrix = {}
        diNucleotides = itertools.product('ACGT', repeat=2)
        for i in diNucleotides:
            matrix[''.join(i)] = []
    elif kind == 'mono':
        matrix = {}
        monoNucleotides = itertools.product('ACGT', repeat=1)
        for i in monoNucleotides:
            matrix[i[0]] = []
    
    lenOfMotif = len(motifs[0])
        
    if kind == 'di':
        for i in matrix.keys():
            matrix[i] = [0]*(lenOfMotif - 1)

        for i in range(lenOfMotif - 1):
            for l in motifs:
                matrix[l[i:i+2]][i] += 1
    elif kind == 'mono':    
        for i in matrix.keys():
            matrix[i] = [0]*lenOfMotif
            
        for i in range(lenOfMotif):
            for l in motifs:
                matrix[l[i]][i] += 1
            
    return(matrix)



def backgroundFreq(seq, kind):
    
    s = ''.join(seq)
    if kind == 'mono':
        background = {}
        monoNucleotides = itertools.product('ACGT', repeat=1)
        for i in monoNucleotides:
            background[i[0]] = s.count(i[0])
    
    elif kind == 'di':
        background = {}
        diNucleotides = itertools.product('ACGT', repeat=2)
        for i in diNucleotides:
            background[''.join(i)] = s.count(i[0])
    
    sumOfNuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sumOfNuc
    
    return(background)


def toScore(normValue, minScore, maxScore):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    score = normValue * (maxScore - minScore) + minScore
    return(score)


def toNorm(score, minScore, maxScore):
    normValue = (score - minScore) / (maxScore - minScore)
    return(normValue)


def minScore(PWM):
    '''
    Вичисляет минимальное значение score для матрицы
    '''
    
    value = int()
    keys = list(PWM.keys())
    lengthPWM = len(PWM[keys[0]])
    for i in range(lengthPWM):
        tmp = []
        for j in keys:
            tmp.append(PWM[j][i])
        value += min(tmp)
    return(value)


def maxScore(PWM):
    '''
    Вичисляет минимальное значение score для матрицы
    '''
    
    value = int()
    keys = list(PWM.keys())
    lengthPWM = len(PWM[keys[0]])
    for i in range(lengthPWM):
        tmp = []
        for j in keys:
            tmp.append(PWM[j][i])
        value += max(tmp)
    return(value)


def score(seq, PWM, kind):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    kind - тип PWM mono or di
    '''
    if kind == 'mono':
        lengthOfSeq = len(seq)
        position = 0
        score = 0
        for letter in seq:
            score += PWM[letter][position]
            position += 1
        return(score)
    
    elif kind == 'di':
        lengthOfSeq = len(seq)
        score = 0
        for i in range(len(seq) - 1):
            twoLetters = seq[i:i+2]
            score += PWM[twoLetters][i]
        return(score)
    else:
        pass


def scanSequences(sequences, PWM, treshold, kind):
    '''
    sequences - список последовательностей (длина последовательностей совпадает с длиной матрицы)
    PWM - матрица
    Сканирование последовательностей матрицей PWM и вычесление score,
    Функция возвращает кол-во последовательностей для которыйх score >= treshold
    kind - тип PWM mono or di
    '''
    aboveTreshold = int()   
    for i in sequences:
        s = score(i, PWM, kind)
        if s >= treshold:
            aboveTreshold += 1
    return(aboveTreshold)


def randomSeq(seq, k):
    '''
    Генерирует k случайных последовательностей с таким же нуклеотидным составом как и seq
    Возвращает list()
    '''
    out = list()
    for i in range(k):
        out.append(''.join(random.sample(seq, k=len(seq))))
    return(out)


def main(fileInput, everyStr, homology, kind):
    '''
    Основная функция
    '''
    
    allMotifs = readFasta(fileInput, everyStr=everyStr)
    allMotifs = removeEqualentSeq(allMotifs, homology=homology)
    results = [] # List([normTreshold_1, FP_1, FN_1],[normTreshold_2, FP_2, FN_2],[] ...) 
    for n in range(len(allMotifs)):
        
        subMotifs = list(allMotifs)
        subMotifs.remove(allMotifs[n])
        
        background = backgroundFreq(subMotifs, kind=kind)
        PCM = makePCM(subMotifs, kind=kind)
        PWM = makePWMfromPCM(PCM, background=background, kind=kind)
        treshold = score(allMotifs[n], PWM, kind=kind)
        
        minS = minScore(PWM)
        maxS = maxScore(PWM)
        normTreshold = toNorm(treshold, minS, maxS)
        
        #Генерация случайных последовательностей
        FP = int()
        totalLen  = int()
        while FP <= 150: 
            randomMotifs = []
            for i in subMotifs:
                randomMotifs += randomSeq(i, 20)
            totalLen += len(randomMotifs) 
            FP += scanSequences(randomMotifs, PWM, treshold, kind=kind)
            if totalLen >= 150000:
                break
        if FP == 0:
            FP = 1.0/(2*totalLen)
        else:
            FP = FP/totalLen
        results.append([normTreshold, FP])
        
    results.sort(key=lambda x: x[1])
    results = results [::-1]
    
    for i in range(len(results)):
        results[i].append((i + 1) / len(allMotifs))
        
    return(results)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='fileInput', required=True , help='Path to file with motifs')
    parser.add_argument('-o', '--output', action='store', dest='fileOutput', required=True, help='Path to file with results of calculation')
    parser.add_argument('-f', '--fasta', action='store_false', dest='fileType', help='Type of input, if file is fasta or file contains [>] use this flag')
    parser.add_argument('-t', '--type', action='store', dest='kind', required=True, help='Type of matrix di- or mononucliotide, after flag print di or mono')
    parser.add_argument('-H', '--removeHomologs', action='store', dest='homology', required=True, help='portion of homology, for example 0.95')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    fileInput = args.fileInput #Путь к фаилу с мотивами
    everyStr = args.fileType #Если строки не разделяются символом '>', то параметр True, иначе пиши False (for Fasta)
    homology = float(args.homology) #Допустимый уровень гомологии между мотивами (если уровень гомологии выше, то последовательность выбрасывается)
    kind = args.kind #Тип матрицы (mono or di)
    fileOutput = args.fileOutput #Путь к фаилу для записи результатов
    
    output = main(fileInput=fileInput, everyStr=everyStr, homology=homology, kind=kind)
    
    with open(fileOutput, 'w') as file:
        file.write('Treshold\tFN rate\tFP rate\n')
        for i in output:
            file.write('{0}\t{2}\t{1}\n'.format(i[0], i[1], i[2]))

