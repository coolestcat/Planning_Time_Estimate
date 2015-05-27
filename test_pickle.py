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

patient_to_timeseries = pickle.load(open("patient_to_timeseries.p", "rb"))

# need to be loaded from nightly
patient_to_priority = pickle.load(open("patient_to_priority.p", "rb"))
patient_to_all_diagnoses = pickle.load(open("patient_to_all_diagnoses.p", "rb"))
patient_to_all_oncologists = pickle.load(open("patient_to_all_oncologists.p", "rb"))
cancer_types_list = pickle.load(open("cancer_types_list.p", "rb"))
onc_list = pickle.load(open("onc_list.p", "rb"))
priorities_list = pickle.load(open("priorities_list.p", "rb"))
qc_holidays = Holidays(country='CA', prov='QC')

# loading the predictors
# svr_pred = pickle.load(open("svr_pred.p", "rb"))
clf_initial = pickle.load(open( "regularized_linear_regression_clf_initial.p", "rb" ))
clf = pickle.load(open( "regularized_linear_regression_clf.p", "rb" ))

print "Done Loading Libraries"

predictions = []
actuals = []
error = False

for i in range(5000):
	patient_id = i
	# for k, v in patient_to_timeseries.items():
		# print k
	v = patient_to_timeseries[patient_id]
	k = patient_id
	for w in v:
		new_feature_vector = []

		if k in patient_to_all_diagnoses.keys():
			for cancer in cancer_types_list:
				if cancer in patient_to_all_diagnoses[k]:
					new_feature_vector.append(float(1.0))
				else:
					new_feature_vector.append(float(0.0))
		else:
			error = True
			break

		if k in patient_to_all_oncologists.keys():
			for oncologist in onc_list:
				if oncologist in patient_to_all_oncologists[k]:
				# if patient_to_oncologist[k]==oncologist:
					new_feature_vector.append(float(1.0))
				else:
					new_feature_vector.append(float(0.0))
		else:
			error = True
			break

		if k in patient_to_priority.keys():
			for priority in priorities_list:
				if patient_to_priority[k]==priority:
					new_feature_vector.append(float(1.0))
				else:
					new_feature_vector.append(float(0.0))
		else:
			error = True
			break

		if k in patient_to_timeseries.keys():
			new_feature_vector.append(float(0.0))
		else: 
			new_feature_vector.append(float(1.0))

		# if w[6] == float(1.0):
		# 	new_feature_vector.append(float(0.0))
		# else:
		# 	new_feature_vector.append(float(1.0))

		# if first==True:
		# 	new_feature_vector.append(float(1.0))
		# else:
		# 	new_feature_vector.append(float(0.0))

		# if (w[5]-w[0]).days < 35: #don't know if final days will be <35 or not! There are some outliers that make the standard deviation much worse
		prediction = clf_initial.predict(new_feature_vector)
		predicted_end_date = w[0] + timedelta(days=prediction)

		others = 0
		real_timedelta = w[5] - w[0]
		half_timedelta = real_timedelta/2 #can you do this?
		half_time = w[0] + half_timedelta

		for a, b in patient_to_timeseries.items():
			for c in b:
				if len(c) > 8:
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
		new_feature_vector.append((w[1]-w[0]).days)
		new_feature_vector.append((w[2]-w[1]).days)
		preds = clf.predict(new_feature_vector)
		# preds = svr_pred.predict(new_feature_vector)
		# print "Prediction: " + str(preds[0])
		# print "Actual: " + str(actualTime(w, qc_holidays))

		predictions.append(preds)
		actuals.append(actualTime(w, 5, 2, qc_holidays))

	if error == True:
		# print "Error - Some Patient Parameter is Missing"
		continue

error_set = []
for i in range(0, len(predictions)):
	error = abs(predictions[i] - actuals[i])
	error_set.append(error)

print "Error: " + str(numpy.mean(error_set))
print "Stddev: " + str(numpy.std(error_set))


# clf = pickle.load(open( "regularized_linear_regression_clf.p", "rb" ))
# testing_set = pickle.load(open("testing_set.p", "rb"))
# target_2 = pickle.load(open("target_2.p", "rb"))

# preds = clf.predict(testing_set)

# error_set = []
# for i in range(0, len(preds)):
# 	error = abs(preds[i] - target_2[i])
# 	error_set.append(error)

# print "Average Absolute Value Error for Regularized Linear Regression: " + str(numpy.mean(error_set))
# print "Stdev of Error: " + str(numpy.std(error_set))

# y_lin = pickle.load(open("svr_y_lin.p", "rb"))

# error_lin = []

# for i in range(0, len(y_lin)):
# 	error2 = abs(y_lin[i] - target_2[i])
# 	# error3 = abs(y_poly[i] - target_2[i])

# 	error_lin.append(error2)
# 	# error_poly.append(error3)

# print "Average Absolute Value Error, linear SVM: " + str(numpy.mean(error_lin))
# print "Stdev of Error, linear SVM: " + str(numpy.std(error_lin))

# svr_days = []
# actual_days = []
# labels = []

# for i in range(50):
# 	# svr_days.append(y_lin[i+50])
# 	# actual_days.append(target_2[i+50])
# 	svr_days.append(y_lin[i+100])
# 	actual_days.append(target_2[i+100])		
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