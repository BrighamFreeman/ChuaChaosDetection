# import drive and needed libraries
# This code was designed to run in Google Colab. You can copy and paste into a Colab file to run
from google.colab import drive
drive.mount("/content/drive", force_remount=True)
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# load in sample csv data sheets
# Data needs to be mounted to drive
shilnikov_data = pd.read_csv("/content/drive/MyDrive/Example_CSV/shilnikov.csv")
shilnikov_data = shilnikov_data.drop_duplicates()
for col in shilnikov_data:
  shilnikov_data[col] = shilnikov_data[col].str.strip()
shilnikov_data = shilnikov_data.rename(columns={"P_1": "alpha", "P_2": "beta", "P_3": 'gamma', 'P_4':
                                                'eigenvalue 1', "P_5": 'eigenvalue 2', 'P_6': 'eigenvalue 3', 'P_7': 'Equilibrium'})
# filter out (0,0,0) solutions
non_zerostring = shilnikov_data[shilnikov_data.Equilibrium != "0+0i"]
#print for debugging
print(len(non_zerostring))

# MATLAB formats complex numbers as a + bi, Numpy formats complex values as a + bj
# Reformat for Python compatibility
shilnikov_data['Equilibrium'] = str(shilnikov_data['Equilibrium'])
shilnikov_data['Equilibrium'] = shilnikov_data['Equilibrium'].str.replace('i','j')
shilnikov_data['Equilibrium'] = np.real(shilnikov_data['Equilibrium'])
print(len(shilnikov_data))

# drop duplicate entries
non_zerostring = non_zerostring.drop_duplicates()

# debugging
non_zerostring['Equilibrium'] = non_zerostring['Equilibrium'].str.replace('i','j').apply(lambda x: np.complex(x))

# Create a new dataframe with rounded values for analysis
rounded_sb = shilnikov_data
rounded_sb['alpha'] = rounded_sb['alpha'].str.replace('i','j').apply(lambda x: complex(x))
rounded_sb['alpha'] = np.real(rounded_sb['alpha'])
rounded_sb['alpha'] = np.around(np.ceil(rounded_sb.alpha.values*2)/2, decimals = 1)


rounded_sb['beta'] = rounded_sb['beta'].str.replace('i','j').apply(lambda x: complex(x))
rounded_sb['beta'] = np.real(rounded_sb['beta'])
rounded_sb['beta'] = np.around(np.ceil(rounded_sb.beta.values*2)/2, decimals = 1)

rounded_sb['gamma'] = rounded_sb['gamma'].str.replace('i','j').apply(lambda x: complex(x))
rounded_sb['gamma'] = np.real(rounded_sb['gamma'])
rounded_sb['gamma'] = np.around(np.ceil(rounded_sb.gamma.values*2)/2, decimals = 1)

# order rounded values by alpha parameter
rounded_sb = rounded_sb.sort_values(by='alpha')

'''

==============================================================
            THIS SECTION IS FOR GRAPHING DATA.
        THIS WAS INTENDED TO RUN IN A JUPYTER NOTEBOOK,
  SO THERE WILL BE COMMENTS DENOTING WHERE CELL BREAKS OCCUR
==============================================================
'''
'''
============================================================
                          NEW CELL
          KDEPlot to graph density of bifurcations
============================================================
'''

sns.set(rc={'figure.figsize':(13,8.27)})
data = [shilnikov_data["alpha"],shilnikov_data["beta"],shilnikov_data["gamma"]]
fig = sns.kdeplot(data=data)
plt.suptitle("Density of Bifurcations for Given Parameter Values")
ax = plt.gca()
ax.set(xlabel="Parameter Value")
plt.show()

'''
================================================================
                          NEW CELL
      Filter KDEPlot to show correlation for alpha and beta
================================================================
'''

data = [shilnikov_data["alpha"],shilnikov_data["beta"]]
fig = sns.kdeplot(data=data)
ax = plt.gca()
ax.set(xlabel="Parameter Value")
plt.suptitle("Concentration of Bifurcations for Beta and Gamma")
plt.show()

'''
================================================================
                          NEW CELL
              Re-round data for easy plotting
================================================================
'''

rounded_sb = shilnikov_data
#rounded_sb['alpha'] = rounded_sb['alpha'].str.replace('i','j').apply(lambda x: np.complex(x))
rounded_sb['alpha'] = np.real(rounded_sb['alpha'])
rounded_sb['alpha'] = np.around(np.ceil(rounded_sb.alpha.values/5)*5, decimals = 0)

#rounded_sb['beta'] = rounded_sb['beta'].str.replace('i','j').apply(lambda x: np.complex(x))
rounded_sb['beta'] = np.real(rounded_sb['beta'])
rounded_sb['beta'] = np.around(np.ceil(rounded_sb.beta.values/5)*5, decimals = 0)

#rounded_sb['gamma'] = rounded_sb['gamma'].str.replace('i','j').apply(lambda x: np.complex(x))
rounded_sb['gamma'] = np.real(rounded_sb['gamma'])
rounded_sb['gamma'] = np.around(np.ceil(rounded_sb.gamma.values/5)*5, decimals = 0)

#display dataframe 
rounded_sb

'''
================================================================
                          NEW CELL
  Get density of bifurcations in different areas of phase space
================================================================
'''

shilnikov_data = shilnikov_data.sort_values(by=['alpha'])
alpha_density = shilnikov_data.value_counts(shilnikov_data['alpha'], ascending=False)
beta_density = shilnikov_data.value_counts(shilnikov_data['beta'])
gamma_density = shilnikov_data.value_counts(shilnikov_data['gamma'])

alpha_df = pd.DataFrame()
alpha_df['count'] = alpha_density
alpha_density

'''
================================================================
                          NEW CELL
Scatterplot to visualize density of bifurcations in phase space
================================================================
'''

sns.scatterplot(alpha_density, marker='s')
sns.scatterplot(beta_density, marker='x')
sns.scatterplot(gamma_density)
plt.suptitle("Density of Bifurcations Around Given Values")
ax = plt.gca()
ax.set(xlabel='Parameter Value', ylabel='Amount of Bifurcations')
plt.legend(labels=['Alpha', 'Beta', 'Gamma'])
alpha_density

'''
================================================================
                          NEW CELL
 Histplot to visualize density of bifurcations in phase space
================================================================
'''

x = len(shilnikov_data['alpha'])
X = np.linspace(-35, 35, num=x)
Y1 = shilnikov_data['alpha']
Y2 = shilnikov_data['beta']
Y3 = shilnikov_data['gamma']

df = pd.DataFrame()
df['alpha'] = Y1
df['beta'] = Y2
df['gamma'] = Y3

ax = plt.gca()
ax.set_xlim([-35.1, 35])
ax.set_ylim([0, 2250])
plt.suptitle("Distribution of Bifurcations Across Parameter Values")
ax.set(xlabel='Parameter Value', ylabel='Amount of Bifurcations')
sns.histplot(df,kde=True, binwidth=5)

'''
================================================================
                          NEW CELL
Regplot to visualize the relationship between variables
================================================================
'''


data_df = pd.DataFrame()
data_df['A_D'] = alpha_density
data_df['B_D'] = beta_density
data_df['C_D'] = gamma_density
data_df['index'] = np.linspace(-150,25, 9)

ind = data_df['index']
Y = data_df['A_D']
Y2 = data_df['B_D']
Y3 = data_df['C_D']

ax = plt.gca()
ax.set_xlim([-150.1, 25.1])
ax.set_ylim([0, 1250])

sns.regplot(data=data_df,x=ind, y=Y, fit_reg=True, label='alpha')
sns.regplot(data=data_df,x=ind, y=Y2, fit_reg=True, label='beta')
sns.regplot(data=data_df,x=ind, y=Y3, fit_reg=True, label='gamma')
plt.legend()
ax.set(xlabel='Parameter Value', ylabel='Amount of Bifurcation Values')
plt.suptitle("Density of Bifurcations Near Parameter Values")

'''
================================================================
                          NEW CELL
        Setting up visualization for Hopf bifurcations
================================================================
'''

#hopf analysis
hopf_data = pd.read_csv("/content/drive/MyDrive/Example_CSV/hopf.csv")
hopf_data = hopf_data.rename(columns={"P_1": "alpha", "P_2": "beta", "P_3": 'gamma', 'P_4': 'eigenvalue 1', "P_5": 'eigenvalue 2', 'P_6': 'eigenvalue 3', 'P_7': 'Previous eigenvalue 1', 'P_8': 'Previous eigenvalue 2', "P_9": 'Previous eigenvalue 3', "P_10": 'Equilibrium', 'P_11': 'Flag'})
hopf_data = hopf_data.drop_duplicates()
hopf_data

eigen_analysis = pd.DataFrame()
eigen_analysis['eigenvalue 1'] = hopf_data['eigenvalue 1']
eigen_analysis['Previous eigenvalue 1'] = hopf_data['Previous eigenvalue 1']
eigen_analysis['eigenvalue 2'] = hopf_data['eigenvalue 2']
eigen_analysis['Previous eigenvalue 2'] = hopf_data['Previous eigenvalue 2']
eigen_analysis['eigenvalue 3'] = hopf_data['eigenvalue 3']
eigen_analysis['Previous eigenvalue 3'] = hopf_data['Previous eigenvalue 3']
eigen_analysis

rounded_hb = hopf_data
rounded_hb['alpha'] = rounded_hb['alpha'].str.replace('i','j').apply(lambda x: complex(x))
rounded_hb['alpha'] = np.real(rounded_hb['alpha'])
rounded_hb['alpha'] = np.around(np.ceil(rounded_hb.alpha.values*2)/2, decimals = 1)
rounded_hb

rounded_hb['beta'] = rounded_hb['beta'].str.replace('i','j').apply(lambda x: complex(x))
rounded_hb['beta'] = np.real(rounded_hb['beta'])
rounded_hb['beta'] = np.around(np.ceil(rounded_hb.beta.values*2)/2, decimals = 1)

rounded_hb['gamma'] = rounded_hb['gamma'].str.replace('i','j').apply(lambda x: complex(x))
rounded_hb['gamma'] = np.real(rounded_hb['gamma'])
rounded_hb['gamma'] = np.around(np.ceil(rounded_hb.gamma.values*2)/2, decimals = 1)

rounded_hb

'''
================================================================
                          NEW CELL
                KDEPlot for Hopf Bifurcations
================================================================
'''

sns.set(rc={'figure.figsize':(13,8.27)})
data = [hopf_data["alpha"],hopf_data["beta"],hopf_data["gamma"]]
fig = sns.kdeplot(data=data)
plt.suptitle("Concentration of Hopf Bifurcations for Given Parameter Values")
ax = plt.gca()
ax.set(xlabel="Parameter Value")
plt.show()

'''
================================================================
                          NEW CELL
                Histplot for Hopf Bifurcations
================================================================
'''


sns.set(rc={'figure.figsize':(13,8.27)})
data = [hopf_data["alpha"],hopf_data["beta"],hopf_data["gamma"]]
fig = sns.histplot(data=data)
plt.suptitle("Concentration of Hopf Bifurcations for Given Parameter Values")
ax = plt.gca()
ax.set(xlabel="Parameter Value")
plt.show()

'''
================================================================
                          NEW CELL
        Set up density dataframe for Hopf bifurcations
================================================================
'''

alpha_density = rounded_hb.value_counts(rounded_hb['alpha'], ascending=False)
beta_density = rounded_hb.value_counts(rounded_hb['beta'])
gamma_density = rounded_hb.value_counts(rounded_hb['gamma'])
data_df = pd.DataFrame()
data_df['A_D'] = alpha_density
data_df['B_D'] = beta_density
data_df['C_D'] = gamma_density
data_df['index'] = np.linspace(-55, 0, 75)
data_df

ind = data_df['index']
Y = data_df['A_D']
Y2 = data_df['B_D']
Y3 = data_df['C_D']

bin1 = data_df['A_D']
bin2 = data_df['B_D'].dropna()
bin3 = data_df['C_D'].dropna()

ax = plt.gca()
ax.set_xlim([-55.2, 0.1])
ax.set_ylim([0, 20])

sns.regplot(data=data_df,x=ind, y=Y, fit_reg=True, scatter_kws = {'s': bin1}, label='alpha')
sns.regplot(data=data_df,x=ind, y=Y2, fit_reg=True, scatter_kws = {'s': bin2}, label='beta')
sns.regplot(data=data_df,x=ind, y=Y3, fit_reg=True, scatter_kws = {'s': bin3}, label='gamma')
ax.set(xlabel='Parameter Value', ylabel='Amount of Bifurcation Values')
plt.legend()
plt.suptitle("Density of Hopf Bifurcations Near Parameter Values")

'''
================================================================
                          NEW CELL
        Import data for runs without gamma parameter
================================================================
'''

shilnikov_noc = pd.read_csv("/content/drive/MyDrive/Example_CSV/shilnikov_nogamma.csv")

shilnikov_noc

shilnikov_noc = shilnikov_noc.rename(columns = {"P_1": 'alpha', "P_2": 'beta', "P_3": 'gamma', "P_4": 'eigenvalue 1', "P_5": 'eigenvalue 2', "P_6": 'eigenvalue 3', "P_7": 'equilibrium'})
shilnikov_noc['alpha'] = shilnikov_noc['alpha'].str.replace('i','j').apply(lambda x: complex(x))
shilnikov_noc['alpha'] = np.real(shilnikov_noc['alpha'])

shilnikov_noc['beta'] = shilnikov_noc['beta'].str.replace('i','j').apply(lambda x: complex(x))
shilnikov_noc['beta'] = np.real(shilnikov_noc['beta'])

shilnikov_noc['gamma'] = shilnikov_noc['gamma'].str.replace('i','j').apply(lambda x: complex(x))
shilnikov_noc['gamma'] = np.real(shilnikov_noc['gamma'])

shilnikov_noc

#add data where gamma = 0 from previous runs

no_gamma = shilnikov_data[shilnikov_data['gamma'] == 0]
shilnikov_noc = shilnikov_noc.append(no_gamma)
shilnikov_noc = shilnikov_noc.drop_duplicates(keep='first')
shilnikov_noc

X = np.linspace(8,20,115)
data = [shilnikov_noc['alpha'], shilnikov_noc['beta']]
print(np.shape(X))
print(np.shape(shilnikov_noc['alpha']))
ax = plt.gca()
sns.kdeplot(data=shilnikov_noc, x=X, y=shilnikov_noc['alpha'], label='alpha')
sns.kdeplot(data=shilnikov_noc, x=X, y=shilnikov_noc['beta'], label='beta')
ax.set(xlabel="beta",ylabel='alpha')
plt.legend()
plt.suptitle("Phase Space Regions with Shilnikov Bifurcation Density")
plt.show()

'''
================================================================
                          NEW CELL
        Import data for runs with cubic non-linearity
================================================================
'''

shilnikov_x3 = pd.read_csv("/content/drive/MyDrive/Example_CSV/shilnikov_x3.csv")
shilnikov_x3 = shilnikov_x3.drop_duplicates()
for col in shilnikov_x3:
  shilnikov_x3[col] = shilnikov_x3[col].str.strip()
shilnikov_x3 = shilnikov_x3.rename(columns={"P_1": "alpha", "P_2": "beta", "P_3": 'gamma', 'P_4': 'eigenvalue 1', "P_5": 'eigenvalue 2', 'P_6': 'eigenvalue 3', 'P_7': 'Equilibrium'})
non_zerostring = shilnikov_x3[shilnikov_x3.Equilibrium != "0+0i"]
print(len(non_zerostring))
shilnikov_x3

shilnikov_x3['alpha'] = shilnikov_x3['alpha'].str.replace('i','j').apply(lambda x: complex(x))
shilnikov_x3['beta'] = shilnikov_x3['beta'].str.replace('i','j').apply(lambda x: complex(x))
shilnikov_x3['gamma'] = shilnikov_x3['gamma'].str.replace('i','j').apply(lambda x: complex(x))

pltdata = [shilnikov_x3['alpha'], shilnikov_x3['beta'], shilnikov_x3['gamma']]
plt.suptitle("KDE Plot of Shilnikov Bifurcations for Cubic Non-Linearity")
sns.kdeplot(data=pltdata)


'''
================================================================
                          NEW CELL
       Rescale data without gamma parameter
================================================================
'''


pltdata = [shilnikov_x3['alpha'], shilnikov_x3['beta']]
plt.suptitle("KDE Plot for Alpha and Beta with Cubic Non-Linearity")
sns.kdeplot(data=pltdata)


'''
================================================================
                          NEW CELL
        Code for alpha and beta parameter histplot
================================================================
'''

ax = plt.gca()
ax.set_xlim([0,1350])
plt.suptitle("Bifurcation Values for alpha and beta")
ax.set(xlabel="Parameter Value")
sns.histplot(data=pltdata, kde=True)

'''
================================================================
                          NEW CELL
        Code for phase space plots of alpha and beta
================================================================
'''
X = np.linspace(8,1500,7600)
data = [shilnikov_x3['alpha'], shilnikov_x3['beta']]
ax = plt.gca()
sns.kdeplot(x=X, y=shilnikov_x3['alpha'], label='alpha')
sns.kdeplot(x=X, y=shilnikov_x3['beta'], label='beta')
ax.set(xlabel="beta",ylabel='alpha')
plt.legend()
plt.suptitle("Phase Space Regions with Shilnikov Bifurcation Density")
plt.show()


