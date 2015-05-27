from sets import Set
from collections import defaultdict
from datetime import datetime
from datetime import date
from datetime import timedelta
from sklearn import decomposition
from sklearn import linear_model
from sklearn import svm
from sklearn import ensemble
from scipy import stats
import numpy.random as nr
import numpy
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

# Patient_ID is the ID in the database
# current stage:
# 0. Ct_Sim but no READY FOR MD CONTOUR
# 1. READY FOR MD CONTOUR but no READY FOR DOSIMETRY
# 2. READY FOR DOSIMETRY but no PRESCRIPTION APPROVED/FAST TRACK 
# 3. PRESCRIPTION APPROVED
# times is an array of all known times. For example, if current_stage is 2, times should be an array of size 3 showing creation dates of:
# [Ct_Sim, READY FOR MD CONTOUR, READY FOR DOSIMETRY]

patient_to_timeseries = pickle.load(open("patient_to_timeseries.p", "rb"))
patient_to_ctstarttimes = pickle.load(open("patient_to_ctstarttimes.p", "rb"))

# need to be loaded from nightly
patient_to_priority = pickle.load(open("patient_to_priority.p", "rb"))
patient_to_all_diagnoses = pickle.load(open("patient_to_all_diagnoses.p", "rb"))
patient_to_all_oncologists = pickle.load(open("patient_to_all_oncologists.p", "rb"))
cancer_types_list = pickle.load(open("cancer_types_list.p", "rb"))
onc_list = pickle.load(open("onc_list.p", "rb"))
priorities_list = pickle.load(open("priorities_list.p", "rb"))
qc_holidays = Holidays(country='CA', prov='QC')

def getEstimate(patient_id, current_stage, times, current_date):
	# all the full timeseries we trained on

	# loading the predictors
	# svr_pred = pickle.load(open("svr_pred.p", "rb"))
	# clf_initial = pickle.load(open( "regularized_linear_regression_clf_initial.p", "rb" ))
	clf = None
	if current_stage == 0:
		clf = pickle.load(open( "regularized_linear_regression_clf.p", "rb" ))
	elif current_stage == 1:
		clf = pickle.load(open( "regularized_linear_regression_clf_1.p", "rb"))
	elif current_stage == 2:
		clf = pickle.load(open( "regularized_linear_regression_clf_2.p", "rb"))
	elif current_stage == 3:
		clf = pickle.load(open( "regularized_linear_regression_clf_3.p", "rb"))
	else:
		# print "Current Stage Parameter is Wrong"
		return -1

	# print "Done Loading Libraries"

	k = patient_id
	new_feature_vector = []

	# Add Diagnoses
	if k in patient_to_all_diagnoses.keys():
		for cancer in cancer_types_list:
			if cancer in patient_to_all_diagnoses[k]:
				new_feature_vector.append(float(1.0))
			else:
				new_feature_vector.append(float(0.0))
	else:
		# print "Patient's Diagnosis not recorded"
		return -1

	# Add Oncologist
	if k in patient_to_all_oncologists.keys():
		for oncologist in onc_list:
			if oncologist in patient_to_all_oncologists[k]:
			# if patient_to_oncologist[k]==oncologist:
				new_feature_vector.append(float(1.0))
			else:
				new_feature_vector.append(float(0.0))
	else:
		# print "Patient's Oncologist not recorded"
		return -1

	# Add Priority
	if k in patient_to_priority.keys():
		for priority in priorities_list:
			if patient_to_priority[k]==priority:
				new_feature_vector.append(float(1.0))
			else:
				new_feature_vector.append(float(0.0))
	else:
		# print "Patient's Priority not recorded"
		return -1

	# If the current treatment course is the first, append 1 else append 0 (First treatment plannings take longer)
	if k in patient_to_timeseries.keys():
		new_feature_vector.append(float(0.0))
	else: 
		new_feature_vector.append(float(1.0))

	# ADD BUSYNESS

	others = 0

	for a, b in patient_to_ctstarttimes.items():
		for c in b:
			if len(c)>1:
				if c[1] == 1:
					if (c[0] < times[0] and c[2] > times[0]) or (c[0] > times[0] and c[0] < current_date):
						others += 1

	# print others
	new_feature_vector.append(others)

	# Add stage times
	for i in range(current_stage):
		new_feature_vector.append((times[i+1]-times[i]).days)

	preds = clf.predict(new_feature_vector)
	# preds = svr_pred.predict(new_feature_vector)[0]
	return preds

	if error == True:
		# print "Error - Some Patient Parameter is Missing"
		return -1

def main():
	# TEST getEstimate
	patient_to_timeseries = pickle.load(open("patient_to_timeseries.p", "rb"))
	qc_holidays = Holidays(country='CA', prov='QC')

	errors = []

	count = 0

	for k, v in patient_to_timeseries.items():
		for this_w in v:

			count += 1
			if count % 500 == 0:
				print count

			times_1 = []
			times_1.append(this_w[0])
			times_1.append(this_w[1])
			times_1.append(this_w[2])
			times_1.append(this_w[3])

			real_timedelta = this_w[5] - this_w[0]
			half_timedelta = real_timedelta/4 #can you do this?
			half_time = this_w[0] + 2*half_timedelta
			cur_date = half_time

			# print "Predicted: " 
			# print getEstimate(k, 3, times_1, cur_date)
			# print "Actual: "
			# print actualTime(this_w, 5, 3, qc_holidays)

			error = abs(actualTime(this_w, 5, 2, qc_holidays) - getEstimate(k, 2, times_1, cur_date))
			if (getEstimate(k, 2, times_1, cur_date) != -1) and (this_w[5]-this_w[0]).days < 35:
				errors.append(error)

	print "Average Absolute Value Error for Regularized Linear Regression: " + str(numpy.mean(errors))
	print "Stdev of Error: " + str(numpy.std(errors))

if __name__ == "__main__":
	main()

