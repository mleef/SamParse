import sys
import math
import locale


locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
WINDOW_SIZE = 250
totalReads = 0
pairedReads = 0
singletonPairs = 0
alignedPairs = 0
notAlignedPairs = 0
sameContPairs = 0
diffContPairs = 0
encounteredReads = {}
contigListSingletons = {}
contigListPairs = {}
diffContigs = {}
contigListSM = {}
pairList = []
insertLengths = []
medians = []
averages = []
standardDevs = []
totalR = 0
totalPR = 0
totalAS = 0
totalNAS = 0
totalAP = 0
totalNAP = 0
totalSCP = 0
totalNCP = 0
pairCoverage = []
singCoverage = []
contigInfo = []
testInfo = []
smCoverage = []
firstFile = True


##Pair to store reads as they are encountered in the file. 
##name: Name of the reads (should be identical)
##contig1: The contig that the first read is on
##contig2: The contig that the second read is on
##position1: Position of first read
##position2: Position of second read
##cigar1: Cigar string of first read
##cigar2: Cigar string of second read
##type1: True if first read has a paired read that aligned, false if it does not.
##type2: True if second read has a paired read that aligned, false if it does not.
class Pair:
	def __init__(self, n, cont1, cont2, pos1, pos2, c1, c2, t1, t2):
		self.name = n
		self.contig1 = cont1
		self.contig2 = cont2
		self.position1 = int(pos1)
		self.position2 = int(pos2)
		self.cigar1 = c1
		self.cigar2 = c2
		self.type1 = t1
		self.type2 = t2
		
	def __str__(self):
		return '({}, {}, {}, {}, {}, {}, {}, {}, {})'.format(self.name, self.contig1, self.contig2, self.position1, self.position2, self.cigar1, self.cigar2, self.type1, self.type2)
		
##Modifies Python's sort to order pairs by location
def pairSort(x,y):
	if(x.position1 > y.position1):
		return 1
	elif(x.position1 == y.position1):
		return 0
	else:
		return -1


##Increments the appropriate contig to keep track of clumps of soft masked reads.
##cont: Contig to increment
##start: Starting position of the soft masking
##length: Length of the soft masked sequence
def incrementSM(cont, start, length):
	global contigListSM

	list = contigListSM[cont]
	
	end = 0
	if(start + length > len(list)):
		end = len(list) - start
	else:
		end = start + length
	
	for i in range(start, end):
		list[i] += 1

##Increments a contig's coverage based on an inputted cigar string. So far only responds to "M" and "S" signals.
##list: Contig list to increment
##start: Starting position of the read.
##cigar: Cigar string of the read.
##name: Name of the contig (for potential soft mask incrementing)	
def incrementContig(list, start, cigar, name):
	a = [-1]*len(cigar)
	b = [-1]*len(cigar)
	start = int(start)
	i = 0
	##print("Before: ")
	##print(list)
	
	for char in cigar:
		if(char.isdigit()):
			a[i] = char
			i = i + 1
		else:
			b[i] = char
			i = i + 1
			
	i = 0
	
	listSize = len(list)
	
	while(i < len(cigar)):
		number = ""
		while(a[i] is not -1):
			number = number + a[i]
			i = i + 1
		command = b[i]
		number = int(number)
		if(command == "M"):
			if(start + number < listSize):
				for j in range(start, start + number):
					list[j] += 1
				start = start + number
			else:
				for j in range(start, listSize):
					list[j] += 1
				break
		elif(command == "S"):
			incrementSM(name, start, number)
			start = start + number
		elif(command == "D"):
			listSize = listSize - number
		elif(command == "I"):
			if(listSize + number > len(list)):
				listSize = len(list)
			else:
				listSize = listSize + number
			start = start + number
		i = i + 1
		
	##print("After: ")
	##print(list)
	
	
		
##Calculates average pair coverage within a window over the length of an entire contig. Stores
##results in "pairCoverage" and the affected contig in "contigInfo"
##name: Name of the contig
##list: List associated with the given contig.
def slidingWindowPairs(name,list):
	##f = open(name, "w")
	##f.write("Region,Coverage")
	##f.write("\n")
	global contigInfo
	global WINDOW_SIZE
	
	
	global pairCoverage
	
	i = 0
	window = WINDOW_SIZE
	curWindow = []
	while(i < len(list)):
		if(i + 500 > len(list)):
			window = len(list) - i
		##print("pairs called with  " + name + " and a window size of " + str(window))
	
		for i in range(i, i + window):
			curWindow.append(list[i])
		##f.write(str(i) + "," + str(sum(curWindow)/window))
		##f.write("\n")
		i = i + 1
		pairCoverage.append(sum(curWindow)/window)
		contigInfo.append(name)

		curWindow[:] = []
	##f.close()
	
	
	
##Calculates average singleton coverage within a window over the length of an entire contig. Stores
##results in "singCoverage".
##name: Name of the contig
##list: List associated with the given contig.		
def slidingWindowSingletons(name,list):
	##f = open(name, "w")
	##f.write("Region,Coverage")
	##f.write("\n")
	global WINDOW_SIZE
	global singCoverage
	global testInfo
	
	i = 0
	window = WINDOW_SIZE
	curWindow = []
	while(i < len(list)):
		if(i + 500 > len(list)):
			window = len(list) - i
		##print("singles called with  " + name + " and a window size of " + str(window))
	
		for i in range(i, i + window):
			curWindow.append(list[i])
		##f.write(str(i) + "," + str(sum(curWindow)/window))
		##f.write("\n")
		i = i + 1
		singCoverage.append(sum(curWindow)/window)
		testInfo.append(name)
		curWindow[:] = []
	##f.close()


##Calculates average soft masked coverage within a window over the length of an entire contig. Stores
##results in "smCoverage".
##name: Name of the contig
##list: List associated with the given contig.		
def slidingWindowSM(name,list):
	##f = open(name, "w")
	##f.write("Region,Coverage")
	##f.write("\n")
	
	global smCoverage
	global WINDOW_SIZE
	
	i = 0
	window = WINDOW_SIZE
	curWindow = []
	while(i < len(list)):
		if(i + 500 > len(list)):
			window = len(list) - i
		##print("singles called with  " + name + " and a window size of " + str(window))
	
		for i in range(i, i + window):
			curWindow.append(list[i])
		##f.write(str(i) + "," + str(sum(curWindow)/window))
		##f.write("\n")
		i = i + 1
		smCoverage.append(sum(curWindow)/window)
		curWindow[:] = []
	##f.close()
				
	
##Calls appropriate methods for contig incrementing based on the alignment characteristics of the pair.
##pair: Pair object to use for incrementing
def calculateCoverage(pair):
	global contigListSingletons
	global contigListPairs
	
	
	if(pair.type1 == True and pair.type2 == True):
		##print("Incrementing " + pair.contig1 + " at position " + str(pair.position1) + " using cigar " + pair.cigar1)
		incrementContig(contigListPairs[pair.contig1], pair.position1, pair.cigar1, pair.contig1)
		##print("Incrementing " + pair.contig2 + " at position " + str(pair.position2) + " using cigar " + pair.cigar2)
		incrementContig(contigListPairs[pair.contig2], pair.position2, pair.cigar2, pair.contig2)
			
	elif(pair.type1 == True and pair.type2 == False):
		##print("Incrementing " + pair.contig1 + " at position " + str(pair.position1) + " using cigar " + pair.cigar1)
		incrementContig(contigListSingletons[pair.contig1], pair.position1, pair.cigar1, pair.contig1)
			
	elif(pair.type1 == False and pair.type2 == True):
		##print("Incrementing " + pair.contig2 + " at position " + str(pair.position2) + " using cigar " + pair.cigar2)
		incrementContig(contigListSingletons[pair.contig2], pair.position2, pair.cigar2, pair.contig2)
		
			
##Increments appropriate counters once both reads of a pair have been encountered (alignedPairs, notAlignedPairs, and singletonPairs)			
def incrementCounters(pair):
	global notAlignedPairs
	global singletonPairs
	global alignedPairs
	
	if(pair.type1 and pair.type2):
		alignedPairs += 1
			
	elif(not pair.type1 and not pair.type2):
		notAlignedPairs +=  1
		
	elif(pair.type1 and not pair.type2 or not pair.type1 and pair.type2):
		singletonPairs +=  1
	
	
				
		


##Writes output statistics as well as contig mappings to separate files. 
def printResults():
	global totalReads
	global pairedReads
	global notAlignedPairs
	global singletonPairs
	global alignedPairs
	global sameContPairs
	global diffContPairs
	global insertLengths
	global encounteredReads
	global totalR
	global totalPR
	global totalAS
	global totalNAS
	global totalAP
	global totalNAP
	global totalSCP
	global totalNCP
	global medians
	global averages
	global standardDevs
	global diffContigs
	global firstFile
	
	totalPairs = totalReads/2
	encounteredReads = {}
	
	outputFile = open("OutputStats", "a")
		

	
	outputFile.write("Total Reads:                " + format(totalReads,"n"))
	outputFile.write("\n")
	outputFile.write("Total Pairs:                " + format(totalPairs,"n"))
	outputFile.write("\n")
	outputFile.write("Aligned Pairs:              " + format(alignedPairs,"n") + "---" + str(float(alignedPairs)/float(totalPairs)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Unaligned Pairs:            " + format(notAlignedPairs,"n") + "---" + str(float(notAlignedPairs)/float(totalPairs)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Singletons:                 " + format(singletonPairs,"n") + "---" + str(float(singletonPairs)/float(totalPairs)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("**********")
	outputFile.write("\n")
	outputFile.write("Pairs on the Same Contig:   " + format(sameContPairs,"n") + "---" +  str(float(sameContPairs)/float(totalPairs)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Pairs on Different Contigs: " + format(diffContPairs,"n") + "---" +  str(float(diffContPairs)/float(totalPairs)*100)[:4] + "%")
	outputFile.write("\n")
	
	count = 0
	f = open("ContigMappings", "a")
	f.write("\n")
	f.write("Name,Frequency")
	f.write("\n")
	for w in sorted(diffContigs, key=diffContigs.get, reverse=True):
		if(count > 99):
			break
		f.write(str(w).replace(",",""))
		f.write(",")
		f.write(str(diffContigs[w]))
  		f.write("\n")
  		count += 1
  		
  	f.close()
  		
	outputFile.write("**********")
	outputFile.write("\n")
	insertLengths.sort()
	length = len(insertLengths)

	if(length%2==0):
		median = insertLengths[length/2]
		medians.append(median)
		outputFile.write("Median Distance Between Reads:  " + str(median)[:4])
		outputFile.write("\n")
	else:
		upper = insertLengths[length/2]
		lower = insertLengths[(length/2)-1]
		median = float(lower + upper)/2
		medians.append(median)
		outputFile.write("Median Distance Between Reads:  " + str(float(lower + upper)/2)[:4])
		outputFile.write("\n")
	
	
	average = sum(insertLengths)/float(len(insertLengths))
	averages.append(average)
	outputFile.write("Average Distance Between Reads: " + str(average)[:4])
	outputFile.write("\n")
	
	tempList = []
	for value in insertLengths:
		tempList.append(math.pow((value - average), 2))
	standardDev = math.sqrt(sum(tempList)/float(len(tempList)))
	standardDevs.append(standardDev)
	outputFile.write("Standard Deviation:             " + str(standardDev)[:4])
	outputFile.write("\n")
	outputFile.close()
	
	
	totalR = totalR + totalReads
	totalPR = totalPR + totalPairs
	totalAS = totalAS + singletonPairs
	totalAP = totalAP + alignedPairs
	totalNAP = totalNAP + notAlignedPairs
	totalSCP = totalSCP + sameContPairs
	totalNCP = totalNCP + diffContPairs
	
	
	totalReads = 0
	pairedReads = 0
	notAlignedPairs = 0
	alignedPairs = 0
	singletonPairs = 0
	sameContPairs = 0
	diffContPairs = 0
	insertLengths[:] = []
	pairList[:] = []
	diffContigs.clear()
	firstFile = False

##Takes in individual lines from the SAM file. Creates Pair objects when reads are initially
##encountered. Once second read is encountered, the appropriate counters are incremented and 
##the coverage is calculated, before the object is discarded. 
##list: Raw line from the SAM file
def processLine(list):
	global insertLengths
	global totalReads
	global encounteredReads
	global pairList
	global pairedReads
	global alignedPairs
	global notAlignedPairs
	global singletonPairs
	global sameContPairs
	global diffContPairs
	global diffContigs
	

	
	name = list[0]
	flag = int(list[1])
	contig = list[2]
	startLoc = list[3]
	cigar = list[5]
	endLoc = list[7]
	
	

	mapped = False
	partnerMapped = False
	properlyAligned = False
		
	if((flag & 4) != 4):
		mapped = True
	if((flag & 8) != 8):
		partnerMapped = True
		
	
	if(mapped and partnerMapped):
		newPair = Pair(name, contig, 0, startLoc, endLoc, cigar, 0, True, True)
		if(name not in encounteredReads):
			encounteredReads[name] = newPair
			##alignedPairs = alignedPairs + 1
		else:
			encounteredReads[name].cigar2 = cigar
			encounteredReads[name].contig2 = contig
			newPair = encounteredReads[name]
			calculateCoverage(newPair)
			incrementCounters(newPair)
			
			if(newPair.contig1 == newPair.contig2):
				sameContPairs += 1
				length = newPair.position2 + len(list[9]) - newPair.position1
				insertLengths.append(length)
			else:
				diffContPairs += 1
				if(newPair.contig1 < newPair.contig2):
					contigPair = (newPair.contig1, newPair.contig2)
				else:
					contigPair = (newPair.contig2, newPair.contig1)

				if(contigPair not in diffContigs):
					diffContigs[contigPair] = 1
				else:
					diffContigs[contigPair] += 1
				
				
			
			del encounteredReads[name]
		return
		
	elif(not mapped and not partnerMapped):
		newPair = Pair(name, contig, 0, startLoc, endLoc, cigar, 0, False, False)
		if(name not in encounteredReads):
			encounteredReads[name] = newPair
			##notAlignedPairs = notAlignedPairs + 1
		else:
			encounteredReads[name].cigar2 = cigar
			encounteredReads[name].contig2 = contig
			newPair = encounteredReads[name]
			calculateCoverage(newPair)
			incrementCounters(newPair)
			del encounteredReads[name]
		return
		
	elif(not mapped and partnerMapped):
		newPair = Pair(name, contig, 0, startLoc, endLoc, cigar, 0, False, True)
		
		if(name not in encounteredReads):
			encounteredReads[name] = newPair
		else:
			encounteredReads[name].cigar2 = cigar
			encounteredReads[name].contig2 = contig
			newPair = encounteredReads[name]
			calculateCoverage(newPair)
			incrementCounters(newPair)
			del encounteredReads[name]
		return
		
	else:
		newPair = Pair(name,contig, 0, startLoc, endLoc, cigar, 0, True, False)
		
		if(name not in encounteredReads):
			encounteredReads[name] = newPair
			##singletonPairs = singletonPairs + 1
		else:
			encounteredReads[name].cigar2 = cigar
			encounteredReads[name].contig2 = contig
			newPair = encounteredReads[name]
			calculateCoverage(newPair)
			incrementCounters(newPair)
			del encounteredReads[name]
		return
		
	


	
			
## For the first file, lists of contigs are read in to be used throughout the program
def parseFile(file):
	global totalReads
	global contigListPairs
	global contigListSingletons
	global contigListSM
	global contSings
	global contPairs
	global firstFile
	
	
	f = open(file)

	while 1:
		line = f.readline()
		if not line: break
		headerLine = line.find("SN:")
		wordList = line.split("\t")
		if(len(wordList) > 2 and headerLine == -1):
			totalReads = totalReads + 1
			processLine(wordList)
		elif(firstFile):
			if(len(wordList) > 1 and headerLine):
				contigName = wordList[1]
				contigLength = wordList[2]
				contigName = contigName[3:]
				contigLength = int(contigLength[3:])
				list = [0 for i in range(contigLength-1)]
				list2 = [0 for i in range(contigLength-1)]
				list3 = [0 for i in range(contigLength-1)]
				if(contigName not in contigListSingletons):
					contigListSingletons[contigName] = list
				if(contigName not in contigListPairs):
					contigListPairs[contigName] = list2
				if(contigName not in contigListSM):
					contigListSM[contigName] = list3
			
				
				
				
	firstFile = False 	
	printResults()
	
	
	##print("# of asterisks in third column: " + str(asteriskCounter))
		


def main(argv = sys.argv):
	global totalR
	global totalPR
	global totalAS
	global totalNAS
	global totalAP
	global totalNAP
	global totalSCP
	global totalNCP
	global medians
	global averages
	global standardDevs
	global pairList
	global contigListSingletons
	global contigListPairs
	global pairCoverage
	global singCoverage
	global smCoverage
	global contigInfo
	global testInfo
	
	outputFile = open("OutputStats", "a")
	f = open("ContigMappings", "a")

	for n in sys.argv[1:]:
		outputFile.write("\n")
		outputFile.write("Input File: " + n)
		outputFile.write("\n")
		f.write("\n")
		f.write("Input File: " + n)
		f.close()
		outputFile.close()
		parseFile(n)
	
	
	for key, value in contigListPairs.iteritems():
		#if(len(value) > 5000):
		slidingWindowPairs(key,value)
		
	
	for key, value in contigListSingletons.iteritems():
		##if(len(value) > 5000):
		slidingWindowSingletons(key,value)
			
	for key, value in contigListSM.iteritems():
		##if(len(value) > 5000):
		slidingWindowSM(key,value)
		
	
			
	allContigs = open("WindowCoverages", "a")
	allContigs.write("Contig, Pair, Singleton, SM, Test")
	allContigs.write("\n")
	
	specialConts = open("ContigsOfInterest", "a")
	seenConts = []
	for i in range(0, len(pairCoverage)):
		if(i > 1 and i < len(pairCoverage)-1):
			if(pairCoverage[i-1] > singCoverage[i-1] and singCoverage[i] > pairCoverage[i] and pairCoverage[i+1] > singCoverage[i+1]):
				if(contigInfo[i] not in seenConts):
					specialConts.write(contigInfo[i])
					specialConts.write("\n")
					seenConts.append(contigInfo[i])
			
		allContigs.write(contigInfo[i] + "," + str(pairCoverage[i]) + "," + str(singCoverage[i]) + "," + str(smCoverage[i]) + "," + testInfo[i])
		allContigs.write("\n")
	
	allContigs.close()
	specialConts.close()
	
	
	outputFile = open("OutputStats", "a")	
	outputFile.write("\n")	
	outputFile.write("TOTALS FOR ALL FILES: ")
	outputFile.write("\n")
	outputFile.write("Total Reads:                      " + format(totalR,"n"))
	outputFile.write("\n")
	outputFile.write("Total Pairs:                      " + format(totalPR,"n"))
	outputFile.write("\n")
	outputFile.write("Total Aligned Pairs:              " + format(totalAP,"n") + "---" + str(float(totalAP)/float(totalPR)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Total Unaligned Pairs:            " + format(totalNAP,"n") + "---" + str(float(totalNAP)/float(totalPR)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Total Aligned Singletons:         " + format(totalAS,"n") + "---" + str(float(totalAS)/float(totalPR)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("**********")
	outputFile.write("\n")
	outputFile.write("Total Pairs on the Same Contig:   " + format(totalSCP,"n") + "---" +  str(float(totalSCP)/float(totalPR)*100)[:4] + "%")
	outputFile.write("\n")
	outputFile.write("Total Pairs on Different Contigs: " + format(totalNCP,"n") + "---" +  str(float(totalNCP)/float(totalPR)*100)[:4] + "%")
	outputFile.write("\n")
	
	median = sum(medians)/float(len(medians))
	average = sum(averages)/float(len(averages))
	standardDev = sum(standardDevs)/float(len(standardDevs))
	outputFile.write("**********")
	outputFile.write("\n")
	outputFile.write("Average Median Value:             " + str(median)[:4])
	outputFile.write("\n")
	outputFile.write("Average Mean Value:               " + str(average)[:4])
	outputFile.write("\n")
	outputFile.write("Average Standard Deviation Value: " + str(standardDev)[:4])
	outputFile.write("\n")
	outputFile.close()
	


	
	

if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass


