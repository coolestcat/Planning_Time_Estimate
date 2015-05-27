import MySQLdb
from sets import Set
from collections import defaultdict
import re
from datetime import datetime
from datetime import date
from datetime import timedelta
from dateutil.easter import easter
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta as rd
from dateutil.relativedelta import MO, TH, FR
from sklearn import decomposition
from sklearn import linear_model
from sklearn import svm
from sklearn import ensemble
from scipy import stats
import numpy.random as nr
import numpy
from sklearn import svm
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import RFECV
import matplotlib.pyplot as plt
from holidays import Holidays
import pickle

def actualTime(w, last, first, qc_holidays):
	# remove weekends
	total_time = (w[last]-w[first]).total_seconds()/float(60*60*24)
	# print "Before: " + str(total_time)
	t_0 = w[first]
	t_5 = w[last]
	to_subtract = 0
	while (t_0 < t_5):
		if (t_0.weekday()==5 or t_0.weekday()==6 or (t_0.date() in qc_holidays)):
			to_subtract += 1
		t_0 += timedelta(days=1)

	total_time = total_time - float(to_subtract)
	if (total_time < 0):
		total_time = 0.0

	return total_time

def printPatient(arr, of, timelines):
	of.write("-----\n")
	of.write("AliasName" + "\t" + "PatientSerNum" + "\t" + "CreationDate" + "\t" + "Status" + "\n")

	for line in arr:
		of.write(str(line["AliasName"]) + "\t" + str(line["PatientSerNum"]) + "\t" + str(line["CreationDate"]) + "\t" + str(line["Status"]))
		for timeline in timelines:
			of.write("\t")
			for dt in timeline:
				if dt == line["CreationDate"]:
					of.write("=====")
		of.write("\n")

	of.write("-----\n\n")

def timelineActive(arr, time):
	toRet = False
	toBreak = False
	for i in range(len(arr)):
		if arr[i]['AliasName'] == "READY FOR MD CONTOUR":
			for j in range(i, len(arr)):
				if arr[i]['AliasName'] == "READY FOR TREATMENT":
					if arr[i]['CreationDate'] < time and arr[j]['CreationDate'] > time:
						toRet = True
						toBreak = True
						break

			if toBreak == True:
				break

	return toRet


def findNext(arr, nextString, iterStart, iterEnd, status):
	l = []
	if status==0:		
		for i in range(iterStart, iterEnd):
			if arr[i]['AliasName'] == nextString:
				l.append(i)
	elif status==1:
		for i in range(iterStart, iterEnd):
			if arr[i]['AliasName'] == nextString and arr[i]['Status']=='Open':
				l.append(i)
	elif status==2:
		for i in range(iterStart, iterEnd):
			if arr[i]['AliasName'] == nextString and arr[i]['Status']=='In Progress':
				l.append(i)
	elif status==3:
		for i in range(iterStart, iterEnd):
			if arr[i]['AliasName'] == nextString and arr[i]['Status']=='Completed':
				l.append(i)		
	elif status==4:
		for i in range(iterStart, iterEnd):
			if arr[i]['AliasName'] == nextString and arr[i]['Status']=='Manually Completed':
				l.append(i)			
	else:
		print "invalid status" # do nothing

	# return list of next positions
	return l

def notTermPresent(arr, notString, start, end, status):
	if status==0:
		for i in range(start, end):
			if arr[i]['AliasName']==notString:
				return 1
	elif status==1:
		for i in range(start, end):
			if arr[i]['AliasName']==notString and arr[i]['Status']=='Open':
				return 1
	elif status==2:
		for i in range(start, end):
			if arr[i]['AliasName']==notString and arr[i]['Status']=='In Progress':
				return 1
	elif status==3:
		for i in range(start, end):
			if arr[i]['AliasName']==notString and arr[i]['Status']=='Completed':
				return 1
	elif status==4:
		for i in range(start, end):
			if arr[i]['AliasName']==notString and arr[i]['Status']=='Manually Completed':
				return 1
	else:
		print "invalid status" # do nothing		

	return 0

def findAll(arr, sequence, not_array, status_array, endCurPatient, l):
	if (len(sequence)==0):
		return l # return [8]

	ret = []

	if not_array[0] == 1: # not simpleExp or functionExp
		newseq =  sequence[:]
		to_check = newseq.pop(0)
		new_not = not_array[:]
		new_not.pop(0)
		new_stat = status_array[:]
		new_stat.pop(0)

		a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)	
		for seq in a:
			if notTermPresent(arr, to_check, seq[0]+1, seq[1], status_array[0])==0:
				ret.append(seq)
	elif not_array[0] == 2: # brace/pipe
		# same if len(l) == 0
		for bracei in range(len(sequence[0])):
			# print braceTerm
			newseq = sequence[:]
			newseq.pop(0)
			newseq.insert(0, sequence[0][bracei])
			new_not = not_array[:]
			new_not.pop(0)
			new_not.insert(0,0)
			new_stat = status_array[:]

			stat_term = new_stat.pop(0)
			new_stat.insert(0, stat_term[bracei])

			a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)
			for seq in a:
				ret.append(seq)	

		# for braceTerm in sequence[0]:
		# 	# print braceTerm
		# 	newseq = sequence[:]
		# 	newseq.pop(0)
		# 	newseq.insert(0, braceTerm)
		# 	new_not = not_array[:]
		# 	new_not.pop(0)
		# 	new_not.insert(0,0)
		# 	new_stat = status_array[:]

		# 	a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)
		# 	for seq in a:
		# 		ret.append(seq)
	elif not_array[0] == 3: # bracestar
		for i in range(len(l)):
			newseq = sequence[:]
			newseq.pop(0)
			new_not = not_array[:]
			new_not.pop(0)
			new_stat = status_array[:]
			new_stat.pop(0)

			a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)
			# print "a: " + str(a)
			for seq in a:
				nextPos = seq[1]
				nar = []
				for bracei in range(len(sequence[0])):
					newl = findNext(arr, sequence[0][bracei], l[i]+1, nextPos, status_array[0][bracei])
					# print "newl: " + str(newl)
					for num in newl:
						nar.append(num)

				# for braceTerm in sequence[0]:
				# 	newl = findNext(arr, braceTerm, l[i], nextPos, status_array[0])
				# 	# print "newl: " + str(newl)
				# 	for num in newl:
				# 		nar.append(num)

				if l[i] < nextPos:
					nar.sort()
					newarr = [l[i]] + nar
					newarr = newarr + seq[1:len(seq)]
					# print "newarr: " + str(newarr)
					ret.append(newarr)
	elif not_array[0] == 4: # star
		if len(l)==0:
			l = []
			for i in range(endCurPatient):
				l.append(i)
			sequence.pop(0)
			not_array.pop(0)
			status_array.pop(0)
			return findAll(arr, sequence, not_array, status_array, endCurPatient, l)

		for i in range(len(l)):
			newseq = sequence[:]
			newseq.pop(0)
			new_not = not_array[:]
			new_not.pop(0)
			new_stat = status_array[:]
			new_stat.pop(0)

			a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)

			for seq in a:
				second = seq[1]
				nar = []

				if l[i] < second:
					for j in range(l[i]+1, second):
						newarr = [l[i]] + [j] + seq[1:len(seq)]
						ret.append(newarr)
	elif not_array[0] == 5: # not brace/pipe/bracestar
		newseq =  sequence[:]
		to_check = newseq.pop(0)
		new_not = not_array[:]
		new_not.pop(0)
		new_stat = status_array[:]
		to_check_stat = new_stat.pop(0)

		a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)	
		for seq in a:
			flag = 0

			for bracei in range(len(to_check)):
				if notTermPresent(arr, to_check[bracei], seq[0]+1, seq[1], to_check_stat[bracei])==1:
					flag = 1
					break

			# for braceTerm in to_check:
			# 	if notTermPresent(arr, braceTerm, seq[0]+1, seq[1], status_array[0])==1:
			# 		flag = 1
			# 		break

			if flag == 0:		
				ret.append(seq)
	elif not_array[0] == 6: # not star
		newseq =  sequence[:]
		to_check = newseq.pop(0)
		new_not = not_array[:]
		new_not.pop(0)
		new_stat = status_array[:]
		new_stat.pop(0)

		a = findAll(arr, newseq, new_not, new_stat, endCurPatient, l)	
		for seq in a:
			if seq[1]==seq[0]+1:
				ret.append(seq)
	else: # normal
		if len(l)==0:
			l = findNext(arr, sequence[0], 0, len(arr), status_array[0])
			sequence.pop(0)
			status_array.pop(0)
			not_array.pop(0)
			if len(l)==0:
				return []
			else:
				return findAll(arr, sequence, not_array, status_array, endCurPatient, l)

		for i in range(len(l)):
			newl = findNext(arr, sequence[0], l[i]+1, endCurPatient, status_array[0]) #[6] #[8]
			# print "newl: " + str(newl)

			if len(newl)!=0:
				newseq = sequence[:]
				newseq.pop(0) # [ready_for_md_contour] []
				new_not = not_array[:]
				new_not.pop(0)
				new_stat = status_array[:]
				new_stat.pop(0)

				a = findAll(arr, newseq, new_not, new_stat, endCurPatient, newl) # a= [6,8]
				# print "a: " + str(a)

				if len(a) > 0:
					if isinstance(a[0], int):
						# print "a: " + str(a)
						for n in a:
							nar = []
							nar.append(l[i])
							nar.append(n)
							ret.append(nar)
					else:
						for seq in a:
							seq.insert(0, l[i]) 
							ret.append(seq)

	# print "ret: " + str(ret)
	return ret 

def main():
	print "Done Loading Libraries"
	mysql_cn = MySQLdb.connect(host='localhost', port=3306, user='alvin', passwd='', db='hig')
	sequenceCursor = mysql_cn.cursor(MySQLdb.cursors.DictCursor)
	
	sequenceCursor.execute("""
    SELECT a.AliasName, t.PatientSerNum, t.CreationDate, t.CompletionDate, t.Status
    FROM Task t INNER JOIN Patient p ON p.PatientSerNum = t.PatientSerNum 
        INNER JOIN Alias a ON a.AliasSerNum = t.AliasSerNum 
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, ap.PatientserNum, ap.ScheduledStartTime, ap.ScheduledEndTime, ap.Status
    FROM Appointment ap INNER JOIN Patient p ON ap.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = ap.AliasSerNum
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, p.PatientSerNum, p.PlanCreationDate, p.PlanCreationDate, p.Status
    FROM Plan p INNER JOIN Patient pa ON p.PatientSerNum = pa.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    WHERE pa.PatientSerNum<35384
    UNION ALL
    SELECT DISTINCT a.AliasName, p.PatientSerNum, t.TreatmentDateTime, p.PlanCreationDate, p.Status 
    FROM TreatmentFieldHstry t INNER JOIN Plan p ON t.PlanSerNum = p.PlanSerNum
    INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    INNER JOIN Patient pat ON pat.PatientSerNum = p.PatientSerNum
    WHERE pat.PatientSerNum<35384
    UNION ALL
    SELECT a.AliasName, d.PatientSerNum, d.ApprovedTimeStamp, d.DateOfService, d.ApprovalStatus
    FROM Document d INNER JOIN Patient p ON d.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a on d.AliasSerNum = a.AliasSerNum
    WHERE p.PatientSerNum<35384
    ORDER BY PatientSerNum, CreationDate
	""")
	
	full_set = []

	# sequence = ['Ct-Sim', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR DOSE CALCULATION', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA' 'READY FOR TREATMENT'], 'READY FOR MD CONTOUR', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR DOSE CALCULATION', ['Ct-Sim', 'Consult', 'Consult Appointment','READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Consult', 'Consult Appointment','Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	# not_array=[0, 5, 0, 5, 0, 5, 2, 5, 0, 5, 0]
	# status=[0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]

	# sequence = ['Ct-Sim', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR DOSE CALCULATION', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA' 'READY FOR TREATMENT'], 'READY FOR MD CONTOUR', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR DOSE CALCULATION', ['Ct-Sim', 'Consult', 'Consult Appointment','READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Consult', 'Consult Appointment','Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	# not_array=[0, 5, 0, 5, 0, 5, 2, 5, 0, 5, 0]
	# status=[0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]

	sequence = ['Ct-Sim', ['Ct-Sim', 'READY FOR DOSE CALCULATION', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA' 'READY FOR TREATMENT'], 'READY FOR MD CONTOUR', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR DOSE CALCULATION', ['Ct-Sim', 'READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	not_array=[0, 5, 0, 5, 0, 5, 2, 5, 0, 5, 0]
	status=[0, [0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]

	cancer_types = Set(['breast','prostate','lung','colon','thyroid', 'bone', 'ovary', 'tongue', 'Hodgkin\'s', 'leukaemia', 'liver', 'skin', 'brain', 'rectum', 'stomach', 'tonsil'])

	t = sequenceCursor.fetchone()
	serNum = t['PatientSerNum']
	curSerNum = serNum

	this_set = []

	while (t is not None):
		arr = []
		while(curSerNum == serNum):
			# print t
			arr.append(t)
			t = sequenceCursor.fetchone()
			if t is None:
				break
			curSerNum = t['PatientSerNum']

		# Do the Sequence Extraction
		new_sequence = sequence[:]
		new_not_array = not_array[:]
		new_status = status[:]


		l = []
		a = findAll(arr, new_sequence, new_not_array, new_status, len(arr), l)

		a_set = Set([])
		for elem in a:
			a_set.add(str(elem))

		a = []
		for elem in a_set:
			a.append(eval(elem))

		# print a
		for timeline in a:
			s = []
			for elem in timeline:
				s.append(arr[elem])

			this_set.append(s)

		if t is None:
			break
		# print "----"

		serNum = t['PatientSerNum']
		curSerNum = serNum

	full_set = full_set + this_set

	# sequence = ['Ct-Sim', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR DOSE CALCULATION', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA' 'READY FOR TREATMENT'], 'READY FOR MD CONTOUR', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR DOSE CALCULATION', ['Ct-Sim', 'Consult', 'Consult Appointment','READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Consult', 'Consult Appointment','Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	# not_array=[0, 5, 0, 5, 0, 5, 2, 5, 0, 5, 0]
	# status=[0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]

	sequence = ['READY FOR MD CONTOUR', '*', 'Ct-Sim', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR DOSE CALCULATION', ['Ct-Sim', 'Consult', 'Consult Appointment','READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Consult', 'Consult Appointment','Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	not_array=[0, 6, 0, 5, 0, 5, 2, 5, 0, 5, 0]
	status=[0, 0, 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]

	# cancer_types = Set(['breast','prostate','lung','colon','thyroid', 'bone', 'ovary', 'tongue', 'Hodgkin\'s', 'leukaemia', 'liver', 'skin', 'brain', 'rectum', 'stomach', 'tonsil'])

	sequenceCursor.execute("""
    SELECT a.AliasName, t.PatientSerNum, t.CreationDate, t.CompletionDate, t.Status
    FROM Task t INNER JOIN Patient p ON p.PatientSerNum = t.PatientSerNum 
        INNER JOIN Alias a ON a.AliasSerNum = t.AliasSerNum 
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, ap.PatientserNum, ap.ScheduledStartTime, ap.ScheduledEndTime, ap.Status
    FROM Appointment ap INNER JOIN Patient p ON ap.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = ap.AliasSerNum
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, p.PatientSerNum, p.PlanCreationDate, p.PlanCreationDate, p.Status
    FROM Plan p INNER JOIN Patient pa ON p.PatientSerNum = pa.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    WHERE pa.PatientSerNum<35384
    UNION ALL
    SELECT DISTINCT a.AliasName, p.PatientSerNum, t.TreatmentDateTime, p.PlanCreationDate, p.Status 
    FROM TreatmentFieldHstry t INNER JOIN Plan p ON t.PlanSerNum = p.PlanSerNum
    INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    INNER JOIN Patient pat ON pat.PatientSerNum = p.PatientSerNum
    WHERE pat.PatientSerNum<35384
    UNION ALL
    SELECT a.AliasName, d.PatientSerNum, d.ApprovedTimeStamp, d.DateOfService, d.ApprovalStatus
    FROM Document d INNER JOIN Patient p ON d.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a on d.AliasSerNum = a.AliasSerNum
    WHERE p.PatientSerNum<35384
    ORDER BY PatientSerNum, CreationDate
	""")

	t = sequenceCursor.fetchone()
	serNum = t['PatientSerNum']
	curSerNum = serNum

	this_set = []

	while (t is not None):
		arr = []
		while(curSerNum == serNum):
			# print t
			arr.append(t)
			t = sequenceCursor.fetchone()
			if t is None:
				break
			curSerNum = t['PatientSerNum']

		# Do the Sequence Extraction
		new_sequence = sequence[:]
		new_not_array = not_array[:]
		new_status = status[:]


		l = []
		a = findAll(arr, new_sequence, new_not_array, new_status, len(arr), l)

		a_set = Set([])
		for elem in a:
			a_set.add(str(elem))

		a = []
		for elem in a_set:
			a.append(eval(elem))

		# print a
		for timeline in a:
			s = []
			for elem in timeline:
				s.append(arr[elem])

			this_set.append(s)

		if t is None:
			break
		# print "----"

		serNum = t['PatientSerNum']
		curSerNum = serNum

	full_set = full_set + this_set
	
	sequence = ['Ct-Sim', ['Ct-Sim', 'READY FOR DOSE CALCULATION', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA' 'READY FOR TREATMENT'], ['Prescription Document (Fast Track)','PRESCRIPTION APPROVED'], ['Ct-Sim', 'Prescription Document (Fast Track)','PRESCRIPTION APPROVED', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR DOSE CALCULATION', 'READY FOR TREATMENT'], 'READY FOR MD CONTOUR', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR MD CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'],  'READY FOR DOSE CALCULATION', ['Ct-Sim', 'READY FOR DOSE CALCULATION', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'READY FOR PHYSICS QA', 'READY FOR TREATMENT'], 'READY FOR PHYSICS QA', ['Ct-Sim', 'Consult', 'Consult Appointment', 'READY FOR PHYSICS QA', 'READY FOR MD CONTOUR', 'READY FOR CONTOUR', 'PRESCRIPTION APPROVED', 'Prescription Document (Fast Track)','READY FOR DOSE CALCULATION'], 'READY FOR TREATMENT']
	not_array=[0, 5, 2, 5, 0, 5,0, 5, 0, 5, 0]
	status=[0, [0,0,0,0,0,0], [0,0],[0,0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0], 0, [0,0,0,0,0,0], 0, [0,0,0,0,0,0,0,0,0], 0]	

	# cancer_types = Set(['breast','prostate','lung','colon','thyroid', 'bone', 'ovary', 'tongue', 'Hodgkin\'s', 'leukaemia', 'liver', 'skin', 'brain', 'rectum', 'stomach', 'tonsil'])

	sequenceCursor.execute("""
    SELECT a.AliasName, t.PatientSerNum, t.CreationDate, t.CompletionDate, t.Status
    FROM Task t INNER JOIN Patient p ON p.PatientSerNum = t.PatientSerNum 
        INNER JOIN Alias a ON a.AliasSerNum = t.AliasSerNum 
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, ap.PatientserNum, ap.ScheduledStartTime, ap.ScheduledEndTime, ap.Status
    FROM Appointment ap INNER JOIN Patient p ON ap.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = ap.AliasSerNum
    WHERE p.PatientSerNum<35384
    UNION ALL 
    SELECT a.AliasName, p.PatientSerNum, p.PlanCreationDate, p.PlanCreationDate, p.Status
    FROM Plan p INNER JOIN Patient pa ON p.PatientSerNum = pa.PatientSerNum
        INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    WHERE pa.PatientSerNum<35384
    UNION ALL
    SELECT DISTINCT a.AliasName, p.PatientSerNum, t.TreatmentDateTime, p.PlanCreationDate, p.Status 
    FROM TreatmentFieldHstry t INNER JOIN Plan p ON t.PlanSerNum = p.PlanSerNum
    INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
    INNER JOIN Patient pat ON pat.PatientSerNum = p.PatientSerNum
    WHERE pat.PatientSerNum<35384
    UNION ALL
    SELECT a.AliasName, d.PatientSerNum, d.ApprovedTimeStamp, d.DateOfService, d.ApprovalStatus
    FROM Document d INNER JOIN Patient p ON d.PatientSerNum = p.PatientSerNum
        INNER JOIN Alias a on d.AliasSerNum = a.AliasSerNum
    WHERE p.PatientSerNum<35384
    ORDER BY PatientSerNum, CreationDate
	""")

	t = sequenceCursor.fetchone()
	serNum = t['PatientSerNum']
	curSerNum = serNum

	this_set = []

	while (t is not None):
		arr = []
		while(curSerNum == serNum):
			# print t
			arr.append(t)
			t = sequenceCursor.fetchone()
			if t is None:
				break
			curSerNum = t['PatientSerNum']

		# Do the Sequence Extraction
		new_sequence = sequence[:]
		new_not_array = not_array[:]
		new_status = status[:]


		l = []
		a = findAll(arr, new_sequence, new_not_array, new_status, len(arr), l)

		a_set = Set([])
		for elem in a:
			a_set.add(str(elem))

		a = []
		for elem in a_set:
			a.append(eval(elem))

		# print a
		for timeline in a:
			s = []
			for elem in timeline:
				s.append(arr[elem])

			this_set.append(s)

		if t is None:
			break
		# print "----"

		serNum = t['PatientSerNum']
		curSerNum = serNum

	full_set = full_set + this_set

	# write timelines to multimap: patient -> list(timelines)
	patient_to_timeseries = defaultdict(list)

	for line in full_set:
		timeline = []
		for item in line:
			timeline.append(item['CreationDate'])

		patient_to_timeseries[line[0]['PatientSerNum']].append(timeline)

	# write tsv of first 50 timelines

	## test patient_to_timeseries

	# count = 0
	# num_patients = 0
	# for k, v in patient_to_timeseries.items():
	# 	num_patients += 1
	# 	for w in v:
	# 		print w
	# 		count += 1
	# 	print "===="

	# print num_patients
	# print count

	# get patient_to_diagnoses

	diagnosisCursor = mysql_cn.cursor(MySQLdb.cursors.DictCursor)

	diagnosisCursor.execute(""" 
	SELECT diag.PatientSerNum, diag.DiagnosisCreationDate, diag.DiagnosisCode, diag.Description
	FROM Patient p INNER JOIN Diagnosis diag ON p.PatientSerNum = diag.PatientSerNum
	ORDER BY PatientSerNum;
	""")

	patient_to_all_diagnoses = defaultdict(list)

	for row in diagnosisCursor:
		splitDesc = re.split(',| ', row['Description'])
		fields = Set(splitDesc)

		for cancer in cancer_types:
			if cancer in fields:
				patient_to_all_diagnoses[row['PatientSerNum']].append(cancer)

	## test patients_to_all_diagnoses

	# for k, v in patient_to_all_diagnoses.items():
	# 	print k
	# 	print "::"
	# 	for w in v:
	# 		print w

	# 	print "====="

	# get patient_to_oncologist

	doctorCursor = mysql_cn.cursor(MySQLdb.cursors.DictCursor)
	
	doctorCursor.execute("""
	SELECT pd.PatientSerNum, pd.DoctorSerNum, pd.OncologistFlag, pd.PrimaryFlag
	FROM PatientDoctor pd INNER JOIN Doctor d ON d.DoctorSerNum = pd.DoctorSerNum
	ORDER BY pd.PatientSerNum;
	""")

	patient_to_oncologist = dict()
	patient_to_all_oncologists = defaultdict(list)
	all_oncologists = Set([])

	for row in doctorCursor:
		if row['OncologistFlag']==1:
			patient_to_all_oncologists[row['PatientSerNum']].append(row['DoctorSerNum'])
			all_oncologists.add(row['DoctorSerNum'])
		# if row['OncologistFlag']==1 and row['PrimaryFlag']==1:
		# 	patient_to_oncologist[row['PatientSerNum']]=row['DoctorSerNum']
		# 	all_oncologists.add(row['DoctorSerNum'])

	# print all_oncologists

	# get patient_to_years

	patientCursor = mysql_cn.cursor(MySQLdb.cursors.DictCursor)

	patientCursor.execute("""
		SELECT * FROM Patient ORDER BY PatientSerNum;
	""")

	patient_to_years = dict()
	for row in patientCursor:
		if (int(row['DateOfBirth'])!=1970): # remove 1970????
			patient_to_years[row['PatientSerNum']]=(2015-int(row['DateOfBirth']))

	## test patient_to_years

	# for k,v in patient_to_years.items():
	# 	print str(k) + "::" + str(v)

	# valid = 0
	# for k, v in patient_to_timeseries.items():
	# 	if k in patient_to_years.keys():
	# 		valid += 1

	# print valid

	# get patient_to_priority

	priorityCursor = mysql_cn.cursor(MySQLdb.cursors.DictCursor)

	priorityCursor.execute("""
		SELECT * FROM Priority;
	""")

	priorities_list = ['SGAS_P1', 'SGAS_P2', 'SGAS_P3', 'SGAS_P4']
	all_priorities = Set([])
	patient_to_priority = dict()

	patient_to_priority_counts = dict()
	# find max number of counts of priority level
	for row in priorityCursor:
		if row['PatientSerNum'] not in patient_to_priority_counts.keys():
			patient_to_priority_counts[row['PatientSerNum']] = [0,0,0,0]
			for i in range(len(priorities_list)):
				if priorities_list[i]==row['PriorityCode']:
					patient_to_priority_counts[row['PatientSerNum']][i] += 1
					break
		else:
			for i in range(len(priorities_list)):
				if priorities_list[i]==row['PriorityCode']:
					patient_to_priority_counts[row['PatientSerNum']][i] += 1
					break

	for k, v in patient_to_priority_counts.items():
		index = 0
		cur_max = 0
		for i in range(len(v)):
			if v[i] > cur_max:
				cur_max = v[i]
				index = i

		patient_to_priority[k] = priorities_list[index]

	# for k, v in patient_to_priority.items():
	# 	print str(k) + " --> " + str(v)


	print "Done Populating Dictionaries"
	# Predictor

	X_matrix = []
	target_set = []
	target_starts = []
	target_ends = []

	# make onc, cancer_type lists
	onc_list = []
	cancer_types_list = []

	for t in all_oncologists:
		onc_list.append(t)

	# print onc_list

	for t in cancer_types:
		cancer_types_list.append(t)

	qc_holidays = Holidays(country='CA', prov='QC')


	# write tsv of all 3 way combinations of priority, doctor, and diagnosis

	# counter = 0
	# for i in range(len(cancer_types_list)):
	# 	for j in range(i, len(onc_list)):
	# 		# for k in range(j, len(priorities_list)):
	# 		print counter
	# 		fileName = "same_profile_" + str(counter) + ".tsv"
	# 		counter += 1
	# 		of = open(fileName, "w")
	# 		cancer_query = cancer_types_list[i]
	# 		doctor_query = onc_list[j]
	# 		# priority_query = priorities_list[k]
	# 		of.write("Cancer Type: " + str(cancer_types_list[i]) + "\n")
	# 		of.write("Doctor: " + str(onc_list[j]) + "\n")

	# 		sequenceCursor.execute("""
	# 	    SELECT a.AliasName, t.PatientSerNum, t.CreationDate, t.CompletionDate, t.Status
	# 	    FROM Task t INNER JOIN Patient p ON p.PatientSerNum = t.PatientSerNum 
	# 	        INNER JOIN Alias a ON a.AliasSerNum = t.AliasSerNum 
	# 	    WHERE p.PatientSerNum<35384
	# 	    UNION ALL 
	# 	    SELECT a.AliasName, ap.PatientserNum, ap.ScheduledStartTime, ap.ScheduledEndTime, ap.Status
	# 	    FROM Appointment ap INNER JOIN Patient p ON ap.PatientSerNum = p.PatientSerNum
	# 	        INNER JOIN Alias a ON a.AliasSerNum = ap.AliasSerNum
	# 	    WHERE p.PatientSerNum<35384
	# 	    UNION ALL 
	# 	    SELECT a.AliasName, p.PatientSerNum, p.PlanCreationDate, p.PlanCreationDate, p.Status
	# 	    FROM Plan p INNER JOIN Patient pa ON p.PatientSerNum = pa.PatientSerNum
	# 	        INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
	# 	    WHERE pa.PatientSerNum<35384
	# 	    UNION ALL
	# 	    SELECT DISTINCT a.AliasName, p.PatientSerNum, t.TreatmentDateTime, p.PlanCreationDate, p.Status 
	# 	    FROM TreatmentFieldHstry t INNER JOIN Plan p ON t.PlanSerNum = p.PlanSerNum
	# 	    INNER JOIN Alias a ON a.AliasSerNum = p.AliasSerNum
	# 	    INNER JOIN Patient pat ON pat.PatientSerNum = p.PatientSerNum
	# 	    WHERE pat.PatientSerNum<35384
	# 	    UNION ALL
	# 	    SELECT a.AliasName, d.PatientSerNum, d.ApprovedTimeStamp, d.DateOfService, d.ApprovalStatus
	# 	    FROM Document d INNER JOIN Patient p ON d.PatientSerNum = p.PatientSerNum
	# 	        INNER JOIN Alias a on d.AliasSerNum = a.AliasSerNum
	# 	    WHERE p.PatientSerNum<35384
	# 	    ORDER BY PatientSerNum, CreationDate
	# 		""")

	# 		t = sequenceCursor.fetchone()
	# 		serNum = t['PatientSerNum']
	# 		curSerNum = serNum

	# 		while (t is not None):
	# 			arr = []
	# 			while(curSerNum == serNum):
	# 				# print t
	# 				arr.append(t)
	# 				t = sequenceCursor.fetchone()
	# 				if t is None:
	# 					break
	# 				curSerNum = t['PatientSerNum']

	# 			this_pat = arr[0]['PatientSerNum']
	# 			if (this_pat in patient_to_all_diagnoses.keys()) and (this_pat in patient_to_all_oncologists.keys()):
	# 				if (cancer_query in patient_to_all_diagnoses[this_pat]) and (doctor_query in patient_to_all_oncologists[this_pat]):
	# 					if (this_pat in patient_to_priority.keys()):
	# 						of.write("Priority: " + str(patient_to_priority[this_pat]) + "\n")
	# 					printPatient(arr, of, patient_to_timeseries[this_pat])

	# 			if t is None:
	# 				break

	# 			serNum = t['PatientSerNum']
	# 			curSerNum = serNum

	# 		of.close()

	print "Done Writing TSVs"

	# do timelines first patient course or not

	for k, v in patient_to_timeseries.items():
		ct_sim_times = []
		for w in v:
			ct_sim_times.append(w[0])

		min_ct_sim_time = min(ct_sim_times)
		for w in v:
			if w[0] == min_ct_sim_time:
				w.append(1.0)
			else:
				w.append(0.0)

	print "Done First or Not"

	# do timeline "busyness"

	# for k, v in patient_to_timeseries.items():
	# 	for w in v:
	# 		others = 0
	# 		for a, b in patient_to_timeseries.items():
	# 			for c in b:
	# 				if (c[5] > w[0] and c[5] < w[5]) or (c[0] > w[0] and c[0] < w[5]) or (c[0] < w[0] and c[5] > w[5]):
	# 					others += 1

	# 		w.append(others)
			# print others


	# print "Done Busyness"


	# create X_matrix and target_set

	for k, v in patient_to_timeseries.items():
		first = True
		for w in v:
			new_feature_vector = []

			if k in patient_to_all_diagnoses.keys():
				for cancer in cancer_types_list:
					if cancer in patient_to_all_diagnoses[k]:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if k in patient_to_all_oncologists.keys():
				for oncologist in onc_list:
					if oncologist in patient_to_all_oncologists[k]:
					# if patient_to_oncologist[k]==oncologist:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			# if k in patient_to_years.keys():
			# 	new_feature_vector.append(float(patient_to_years[k]))
			# else:
			# 	break

			if k in patient_to_priority.keys():
				for priority in priorities_list:
					if patient_to_priority[k]==priority:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if w[6] == float(1.0):
				new_feature_vector.append(float(0.0))
				# print "1.0"
			else:
				new_feature_vector.append(float(1.0))
				# print "0.0"

			# new_feature_vector.append(float(w[7]))

			# if first==True:
			# 	new_feature_vector.append(float(1.0))
			# else:
			# 	new_feature_vector.append(float(0.0))

			if (w[5]-w[0]).days < 35: #only take one from each patient???
				X_matrix.append(new_feature_vector)
				# remove weekends
				total_time = (w[5]-w[0]).total_seconds()/float(60*60*24)
				# print "Before: " + str(total_time)
				t_0 = w[0]
				t_5 = w[5]
				to_subtract = 0
				while (t_0 < t_5):
					if (t_0.weekday()==5 or t_0.weekday()==6 or (t_0.date() in qc_holidays)):
						to_subtract += 1
					t_0 += timedelta(days=1)

				total_time = total_time - float(to_subtract)
				if (total_time < 0):
					total_time = 0.0

				# print "After: " + str(total_time)

				target_set.append(total_time)

			first = False

	print "Done Building Matrix"
	training_set = [X_matrix[i] for i in range(0, len(X_matrix) - len(X_matrix)/5)]
	testing_set = [X_matrix[i] for i in range(len(X_matrix)-len(X_matrix)/5, len(X_matrix))]

	target_1 = [target_set[i] for i in range(0, len(X_matrix) - len(X_matrix)/5)]
	target_2 = [target_set[i] for i in range(len(X_matrix)-len(X_matrix)/5, len(X_matrix))]
	# for i in range(0,20):
	# 	print str(X_matrix[i]) + " -> " + str(target_set[i])

	print "Number of training samples: " + str(len(X_matrix))
	print "Number of target values: " + str(len(target_set))

	clf = linear_model.RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1000.0])
	clf.fit(X_matrix, target_set)

	pickle.dump(clf, open( "regularized_linear_regression_clf_initial.p", "wb" ))

	# create X_matrix and target_set again
	print "Creating X_Matrix and Target set again // Making All Predictions"
	X_matrix = []
	target_set = []

	for k, v in patient_to_timeseries.items():
		first = True
		for w in v:
			new_feature_vector = []

			if k in patient_to_all_diagnoses.keys():
				for cancer in cancer_types_list:
					if cancer in patient_to_all_diagnoses[k]:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if k in patient_to_all_oncologists.keys():
				for oncologist in onc_list:
					if oncologist in patient_to_all_oncologists[k]:
					# if patient_to_oncologist[k]==oncologist:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			# if k in patient_to_years.keys():
			# 	new_feature_vector.append(float(patient_to_years[k]))
			# else:
			# 	break

			if k in patient_to_priority.keys():
				for priority in priorities_list:
					if patient_to_priority[k]==priority:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if w[6] == float(1.0):
				new_feature_vector.append(float(0.0))
				# print "1.0"
			else:
				new_feature_vector.append(float(1.0))
				# print "0.0"

			# new_feature_vector.append(float(w[7]))

			# if first==True:
			# 	new_feature_vector.append(float(1.0))
			# else:
			# 	new_feature_vector.append(float(0.0))

			if (w[5]-w[0]).days < 35: #only take one from each patient???
				prediction = clf.predict(new_feature_vector)
				predicted_end_date = w[0] + timedelta(days=prediction)
				w.append(predicted_end_date)
				# print len(w)

	print "Finally Adding Business Parameter/Stages:"
	X_matrix = []
	target_set = []

	X_matrix_1 = []
	X_matrix_2 = []
	X_matrix_3 = []
	X_matrix_4 = []
	X_matrix_5 = []

	target_set_1 = []
	target_set_2 = []
	target_set_3 = []
	target_set_4 = []
	target_set_5 = []

	for k, v in patient_to_timeseries.items():
		first = True
		for w in v:
			new_feature_vector = []

			if k in patient_to_all_diagnoses.keys():
				for cancer in cancer_types_list:
					if cancer in patient_to_all_diagnoses[k]:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if k in patient_to_all_oncologists.keys():
				for oncologist in onc_list:
					if oncologist in patient_to_all_oncologists[k]:
					# if patient_to_oncologist[k]==oncologist:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			# if k in patient_to_years.keys():
			# 	new_feature_vector.append(float(patient_to_years[k]))
			# else:
			# 	break

			if k in patient_to_priority.keys():
				for priority in priorities_list:
					if patient_to_priority[k]==priority:
						new_feature_vector.append(float(1.0))
					else:
						new_feature_vector.append(float(0.0))
			else:
				break

			if w[6] == float(1.0):
				new_feature_vector.append(float(0.0))
				# print "1.0"
			else:
				new_feature_vector.append(float(1.0))
				# print "0.0"

			# new_feature_vector.append(float(w[7]))

			# if first==True:
			# 	new_feature_vector.append(float(1.0))
			# else:
			# 	new_feature_vector.append(float(0.0))

			if (w[5]-w[0]).days < 35: #only take one from each patient???
				prediction = clf.predict(new_feature_vector)
				predicted_end_date = w[0] + timedelta(days=prediction)

				others = 0
				real_timedelta = w[5] - w[0]
				half_timedelta = real_timedelta/4 #can you do this?
				half_time = w[0] + 4*half_timedelta

				for a, b in patient_to_timeseries.items():
					for c in b:
						if len(c) > 7:
							if c[5] > half_time:
								if (c[0] < w[0] and c[7] > w[0]) or (c[0] > w[0] and c[0] < half_time):
									others += 1
							else:
								if (c[0] < w[0] and c[5] > w[0]) or (c[0] > w[0] and c[0] < half_time):
									others += 1

						# if (c[5] > w[0] and c[5] < w[5]) or (c[0] > w[0] and c[0] < w[5]) or (c[0] < w[0] and c[5] > w[5]):
						# 	others += 1

				# print others
				w.append(others)
				new_feature_vector.append(others)
				X_matrix.append(new_feature_vector)

				# Ct-Sim Segment Included
				new_feature_vector = new_feature_vector[:]
				new_feature_vector.append((w[1]-w[0]).days)
				X_matrix_1.append(new_feature_vector)

				# MD Contour Segment Included
				new_feature_vector = new_feature_vector[:]
				new_feature_vector.append((w[2]-w[1]).days)
				X_matrix_2.append(new_feature_vector)

				# Dosimetry Segment Included
				new_feature_vector = new_feature_vector[:]
				new_feature_vector.append((w[3]-w[2]).days)
				X_matrix_3.append(new_feature_vector)

				# Prescription Approved Segment Included
				new_feature_vector = new_feature_vector[:]
				new_feature_vector.append((w[4]-w[3]).days)
				X_matrix_4.append(new_feature_vector)

				# Ready for Physics QA Segment Included
				new_feature_vector = new_feature_vector[:]
				new_feature_vector.append((w[5]-w[4]).days)
				X_matrix_5.append(new_feature_vector)

				target_set.append(actualTime(w, 5, 0, qc_holidays))
				target_set_1.append(actualTime(w, 5, 1, qc_holidays))
				target_set_2.append(actualTime(w, 5, 2, qc_holidays))
				target_set_3.append(actualTime(w, 5, 3, qc_holidays))
				target_set_4.append(actualTime(w, 5, 4, qc_holidays))
				# target_set_5.append(actualTime(w, 5, 5, qc_holidays))


			first = False

	X_matrix = X_matrix_2
	target_set = target_set_2


	training_set = [X_matrix[i] for i in range(0, len(X_matrix) - len(X_matrix)/5)]
	testing_set = [X_matrix[i] for i in range(len(X_matrix)-len(X_matrix)/5, len(X_matrix))]

	target_1 = [target_set[i] for i in range(0, len(X_matrix) - len(X_matrix)/5)]
	target_2 = [target_set[i] for i in range(len(X_matrix)-len(X_matrix)/5, len(X_matrix))]
	# for i in range(0,20):
	# 	print str(X_matrix[i]) + " -> " + str(target_set[i])

	print "Number of training samples: " + str(len(X_matrix))
	print "Number of target values: " + str(len(target_set))
	print "Length of feature vector: " + str(len(X_matrix[0]))

	print "ACTUAL PREDICTIONS (ROUND 2):"

	pickle.dump(patient_to_timeseries, open("patient_to_timeseries.p", "wb"))
	pickle.dump(patient_to_priority, open("patient_to_priority.p", "wb"))
	pickle.dump(patient_to_all_diagnoses, open("patient_to_all_diagnoses.p", "wb"))
	pickle.dump(patient_to_all_oncologists, open("patient_to_all_oncologists.p", "wb"))
	pickle.dump(cancer_types_list, open("cancer_types_list.p", "wb"))
	pickle.dump(onc_list, open("onc_list.p", "wb"))
	pickle.dump(priorities_list, open("priorities_list.p", "wb"))

	# print "Ordinary Least Squares Whole Dataset: "
	# clf = linear_model.LinearRegression(normalize = True)

	# clf.fit(X_matrix, target_set)
	# preds = clf.predict(X_matrix)

	# error_set = []

	# for i in range(0, len(X_matrix)):
	# 	error = abs(preds[i] - target_set[i])
	# 	error_set.append(error)

	# average_error = numpy.mean(error_set)
	# stdev = numpy.std(error_set)

	# print "Average Absolute Value Error: " + str(average_error)
	# print "Stdev of Error: " + str(stdev)

	# clf.fit(training_set, target_1)
	# preds = clf.predict(testing_set)

	# total_error = 0.0
	# error_set = []

	# for i in range(0, len(preds)):
	# 	error = abs(preds[i] - target_2[i])
	# 	error_set.append(error)

	# print "Average Absolute Value Error for Ordinary Linear Regression: " + str(numpy.mean(error_set))
	# print "Stdev of Error: " + str(numpy.std(error_set))

	print "Regularized Least Squares: "

	clf = linear_model.RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1000.0])
	clf.fit(training_set, target_1)

	print "Regularization Parameter: " + str(clf.alpha_)

	preds = clf.predict(testing_set)

	error_set = []
	for i in range(0, len(preds)):
		error = abs(preds[i] - target_2[i])
		error_set.append(error)


	print "Average Absolute Value Error for Regularized Linear Regression: " + str(numpy.mean(error_set))
	print "Stdev of Error: " + str(numpy.std(error_set))

	pickle.dump(clf, open( "regularized_linear_regression_clf.p", "wb" ))
	pickle.dump(testing_set, open("testing_set.p", "wb"))
	pickle.dump(target_2, open("target_2.p", "wb"))

	print "Done Pickling"

	# params = {'n_estimators': 700, 'max_depth': 5, 'min_samples_split': 1,
 #          'learning_rate': 0.01, 'loss': 'ls'}
	# clf = ensemble.GradientBoostingRegressor(**params)

	# clf.fit(training_set, target_1)
	# preds = clf.predict(testing_set)

	# total_error = 0.0
	# error_set = []

	# for i in range(0, len(preds)):
	# 	error = abs(preds[i] - target_2[i])
	# 	error_set.append(error)

	# print "Average Absolute Value Error for Gradient Tree Boosting Regressor: " + str(numpy.mean(error_set))
	# print "Stdev of Error: " + str(numpy.std(error_set))

	svr_days = []
	actual_days = []
	labels = []

	for i in range(50):
		# svr_days.append(y_lin[i+50])
		# actual_days.append(target_2[i+50])
		svr_days.append(preds[i+50])
		actual_days.append(target_2[i+50])		
		labels.append(str(i))

	# print "here"
	N = 50
	ind = numpy.arange(N)
	width = 0.35

	fig, ax = plt.subplots()

	rects0 = ax.bar(ind, svr_days, width, color='r')
	rects1 = ax.bar(ind+width, actual_days, width, color='g')

	x = range(50)
	y = []
	for i in x:
		y.append(14)

	ax.plot(x,y)

	ax.set_ylabel('Days')
	ax.set_title('A comparison of our estimate and the actual wait time for the first 50 test cases')
	ax.set_xticks(ind + width)
	ax.set_xticklabels(labels)
	plt.yticks(numpy.arange(0, 22, 1))
	plt.ylim(ymin=0, ymax=22)

	plt.ylim(ymin=0)

	plt.show()

	# pca = decomposition.PCA(n_components=3)
	# pca.fit(training_set)
	# training_set = pca.transform(training_set)
	# testing_set = pca.transform(testing_set)

	# print "Done decomposition with PCA"

	# clf = linear_model.LinearRegression(normalize = True)
	# clf.fit(training_set, target_1)

	# selector = RFECV(clf, step=1, cv=5)
	# selector = selector.fit(training_set, target_1)

	# print "Support and Ranking: " 
	# print selector.support_
	# print selector.ranking_

	# preds = selector.predict(testing_set)

	# error_set = []

	# for i in range(0, len(preds)):
	# 	error = abs(preds[i] - target_2[i])
	# 	error_set.append(error)

	# average_error = numpy.mean(error_set)
	# stdev = numpy.std(error_set)
	# print "Average Absolute Value Error of Ordinary Linear Regression, With Recursive Feature Elimination: " + str(average_error)
	# print "Stdev of Error: " + str(stdev)

	# # Random Forests

	# rand_forest = RandomForestRegressor()

	# y_rand_forest = rand_forest.fit(training_set, target_1).predict(testing_set)


	# errors = []
	# # error_poly = []

	# for i in range(0, len(y_rand_forest)):
	# 	error = abs(y_rand_forest[i] - target_2[i])
	# 	errors.append(error)

	# print "Average Absolute Value Error, Random Forest: " + str(numpy.mean(errors))
	# print "Stdev of Error, Random Forest: " + str(numpy.std(errors))


	# svm

	# svr_lin = svm.SVR(kernel='linear', C=1e3)
	# # svr_poly = svm.SVR(kernel='poly', C=1e3, degree=2)

	# svr_pred = svr_lin.fit(training_set, target_1)
	# y_lin = svr_pred.predict(testing_set)
	# # y_lin = svr_lin.fit(training_set, target_1).predict(testing_set)
	# # y_poly = svr_poly.fit(training_set, target_1).predict(testing_set)

	# error_lin = []
	# # error_poly = []

	# for i in range(0, len(y_lin)):
	# 	error2 = abs(y_lin[i] - target_2[i])
	# 	# error3 = abs(y_poly[i] - target_2[i])

	# 	error_lin.append(error2)
	# 	# error_poly.append(error3)

	# print "Average Absolute Value Error, linear SVM: " + str(numpy.mean(error_lin))
	# print "Stdev of Error, linear SVM: " + str(numpy.std(error_lin))

	# pickle.dump(svr_pred, open( "svr_pred.p", "wb" ))

	# svr_days = []
	# actual_days = []
	# labels = []

	# for i in range(50):
	# 	# svr_days.append(y_lin[i+50])
	# 	# actual_days.append(target_2[i+50])
	# 	svr_days.append(y_lin[i+50])
	# 	actual_days.append(target_2[i+50])		
	# 	labels.append(str(i))

	# N = 50
	# ind = numpy.arange(N)
	# width = 0.35

	# fig, ax = plt.subplots()

	# rects0 = ax.bar(ind, svr_days, width, color='r')
	# rects1 = ax.bar(ind+width, actual_days, width, color='g')

	# x = range(50)
	# y = []
	# for i in x:
	# 	y.append(14)

	# ax.plot(x,y)

	# ax.set_ylabel('Days')
	# ax.set_title('A comparison of our estimate and the actual wait time for the first 50 test cases')
	# ax.set_xticks(ind + width)
	# ax.set_xticklabels(labels)
	# plt.yticks(numpy.arange(0, 22, 1))
	# plt.ylim(ymin=0, ymax=22)

	# plt.ylim(ymin=0)

	# plt.show()


if __name__ == "__main__":
	main()