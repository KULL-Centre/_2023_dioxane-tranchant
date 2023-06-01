import numpy as np
import matplotlib.pyplot as plt
import csv

def chi2_individualerror(dataset):
    '''Calculate the chi2 vs dioxane Rh for a given dataset in the format: D/D, SE(D/D), Rh, SE(Rh).'''
    datalist = []
    with open(dataset, 'r') as datafile:
        data = csv.reader(datafile, delimiter=';')
        for row in data:
            datalist.append(row)

    #Loop that iterates chi2 calculation for each line in the dataset (excluding the header)
    chi2_list = []
    for line in datalist[1:]:
        #Temporary list for storing chi2 values for each protein. Values will be summed after the loop
        chi2_list_temp = []
        for i in np.arange(2.0,2.5,0.01):
            chi2_val = (float(line[1])-(float(line[3])/i))**2/(float(line[2])**2+float(line[4])**2)
            chi2_list_temp.append(chi2_val)
        chi2_list.append(chi2_list_temp)
    chi2_array = np.array(chi2_list)
    chi2_sum = np.sum(chi2_array, 0)

    Rh_diox = np.arange(2.0,2.5,0.01)
    #Plot the summed chi2 value against each assumed dioxane Rh value
    plt.plot(Rh_diox, chi2_sum)

    #Plot each protein set chi2 value individually against each assumed dioxane Rh value
    #plt.plot(Rh_diox, chi2_array[0], label='HEWL'), plt.plot(Rh_diox, chi2_array[1], label='RNaseA'), plt.plot(Rh_diox, chi2_array[2], label='Myoglobin'), plt.plot(Rh_diox, chi2_array[3], label='ACBP'), plt.plot(Rh_diox, chi2_array[4], label='S100A13'), plt.plot(Rh_diox, chi2_array[5], label='Prolactin'), plt.plot(Rh_diox, chi2_array[6], label='14-3-3')

    ymin = min(chi2_sum)
    xpos = np.where(chi2_sum == ymin)
    xmin = Rh_diox[xpos]
    print('The chi2-minimum is at the assumed dioxane Rh of', xmin)

    plt.xlabel(r'$R_h$ of 1,4-dioxane [$\AA$]')
    plt.ylabel(r'$\chi^2$')
    plt.ylim(0, 500)

    plt.show()

chi2_individualerror('chi2_data.csv')
