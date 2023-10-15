import requests
import pandas as pd
import json
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sodapy import Socrata
# Get information about sites

class my_dictionary(dict):
 
  # __init__ function
  def __init__(self):
    self = dict()
 
  # Function to add key:value
  def add(self, key, value):
    self[key] = value

state_names = []
state_ids = ["KS","MA", "TX", "UT","FL", "OK", "ND", "AZ", "WI", "WY", "GA"]


final_datasources = my_dictionary()

for state in state_names:
    STATE = state

    # CDC Queries
    salmon = ["52cr-rw4k",	"d6kj-devz",	"4qb4-rsd8",	"hwyq-75wu",	"6rpz-c2y5",	"fuzh-wm4c",	"mvaf-qxac",	"2vtj-68zm",	"rcdh-n3ej"]
    client = Socrata("data.cdc.gov", None)

    bins = []
    for i in range(1, 14):
        bins = bins + [i] * 4

    bins.append(13)

    for code in salmon:
        results = client.get(code, where = "reporting_area='" + STATE + "'")
        results_df = pd.DataFrame.from_records(results)
        final_datasources.add(code, results_df)

    i = 0

    for key in final_datasources:
        columns_df = list(final_datasources[key].columns)

        selection_base = ['mmwr_year', 'mmwr_week']

        if "salmonella_paratyphi_infection_1" in columns_df:
            if "salmonellosis_excluding" not in columns_df:
                continue
            
            addition_select = ["salmonellosis_excluding"]

        elif "salmonellosis_excluding_paratyphoid_fever_andtyphoid_fever_previous_52_weeks_med" in columns_df:
            if "salmonellosis_excluding_paratyphoid_fever_andtyphoid_fever_current_week" not in columns_df:
                continue
            addition_select = ["salmonellosis_excluding_paratyphoid_fever_andtyphoid_fever_current_week"]

        elif "salmonellosis_excluding_salmonella_paratyphi_infection_and_salmonella_typhoid_infection_current_week" in columns_df:
            if "salmonellosis_excluding_salmonella_paratyphi_infection_and_salmonella_typhoid_infection_current_week" not in columns_df:
                continue
            addition_select = ["salmonellosis_excluding_salmonella_paratyphi_infection_and_salmonella_typhoid_infection_current_week"]

        elif "salmonellosis_current_week" in columns_df:
            if "salmonellosis_current_week" not in columns_df:
                continue

            addition_select = ["salmonellosis_current_week"]
        else:
            continue

        selection_final = selection_base + addition_select

        mini_df = final_datasources[key][selection_final]
        mini_df = mini_df.rename(columns = {addition_select[0]:"salmonellosis"})

        mini_df["mmwr_year"]= pd.to_numeric(mini_df["mmwr_year"])
        mini_df["mmwr_week"]= pd.to_numeric(mini_df["mmwr_week"])
        mini_df["salmonellosis"]= pd.to_numeric(mini_df["salmonellosis"])
        mini_df = mini_df.fillna(0)
        mini_df = mini_df.sort_values("mmwr_week")
        mini_df = mini_df.groupby(["mmwr_year", "mmwr_week"]).agg(sum)
        mini_df = mini_df.reset_index()

        mini_df["bin"] = bins[0:(len(mini_df))]

        mini_df = mini_df.sort_values(by = "mmwr_week")

        if i == 0:
            plotting_df = mini_df
            
            i += 1

        else:
            plotting_df = pd.concat([plotting_df, mini_df])

    plotting_df["mmwr_year"]= pd.to_numeric(plotting_df["mmwr_year"])
    plotting_df["mmwr_week"]= pd.to_numeric(plotting_df["mmwr_week"])
    plotting_df["salmonellosis"]= pd.to_numeric(plotting_df["salmonellosis"])

    plotting_df = plotting_df.fillna(0)


    plotting_df["COVID1"] = plotting_df["mmwr_year"] == 2019
    plotting_df["COVID2"] = plotting_df["mmwr_year"] == 2020
    plotting_df["COVID3"] = plotting_df["mmwr_year"] >= 2021


    X_data = plotting_df[["mmwr_week", "COVID1", "COVID2", "COVID3"]]
    Y_data = plotting_df["salmonellosis"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)

    binned = plotting_df.groupby(["mmwr_year","bin"]).agg(sum)
    binned = binned.reset_index()

    binned["COVID1"] = (binned["mmwr_year"] >= 2020) * 1

    X_data = binned[["bin", "COVID1"]]
    Y_data = binned["salmonellosis"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)

    modifier = {}

    for item in np.unique(bins):
        mini_binned = binned[binned["bin"] == item]
        X_data = mini_binned[["COVID1"]]
        Y_data = mini_binned["salmonellosis"]

        regressor = LinearRegression()

        regressor.fit(X_data, Y_data)
        # Now add coefficients where needed 
        modifier.update({item:regressor.coef_})

    adjustment_column = []
    for row in binned.iterrows():
        current_bin = list(row[1])[1]

        if row[1]["COVID1"] == 1:
            mod_value = list(modifier[current_bin])[0]
        
        else:
            mod_value = 0

        adjustment_column.append(mod_value)

    binned["modifier"] = adjustment_column
    binned["modifier"] = binned["modifier"] * -1
    binned["covid_adjusted"] = binned["salmonellosis"] + binned["modifier"]

    average_plot = binned.groupby("bin").agg(np.mean).reset_index()

    binned['scaled_v3'] = ((binned['covid_adjusted'] - binned.groupby('bin')['covid_adjusted'].transform('mean')) / \
                    (binned.groupby('bin')['covid_adjusted'].transform('std')))
    
    binned.to_csv(STATE+"_CDC_sal.csv")

    codes = ["23gt-ssfe", "ydsy-yh5w", "yqwx-bvu7", "n24i-76tn", "rwap-xbst", "r9mz-pvtk", "mrip-2k2a", "xvdv-hq7x"]

    final_datasources_ip = my_dictionary()

    client = Socrata("data.cdc.gov", None)

    for code in codes:
        results = client.get(code, where = "reporting_area='" + STATE + "'")
        results_df = pd.DataFrame.from_records(results)
        final_datasources_ip.add(code, results_df)

    bins = []
    for i in range(1, 14):
        bins = bins + [i] * 4

    bins.append(13)

    i = 0
    for key in final_datasources_ip:
        columns_df = list(final_datasources_ip[key].columns)

        selection_base = ['mmwr_year', 'mmwr_week']

        if "invasive_pneumococcal_disease_all_ages_current_week" in columns_df:

            additional = ["invasive_pneumococcal_disease_all_ages_current_week"]

        elif "invasive_pneumococcal_disease_all_ages_confirmed_current_week" in columns_df:
            additional = ["invasive_pneumococcal_disease_all_ages_confirmed_current_week"]
        
        elif "invasive_pneumococcal_disease" in columns_df:
            additional = ["invasive_pneumococcal_disease"]

        else:
            continue

        selection_final = selection_base + additional

        mini_df = final_datasources_ip[key][selection_final]
        mini_df = mini_df.rename(columns = {additional[0]:"invasive_pneumococcal"})
        
        mini_df["mmwr_year"]= pd.to_numeric(mini_df["mmwr_year"])
        mini_df["mmwr_week"]= pd.to_numeric(mini_df["mmwr_week"])
        mini_df["invasive_pneumococcal"]= pd.to_numeric(mini_df["invasive_pneumococcal"])
        mini_df = mini_df.fillna(0)
        mini_df = mini_df.sort_values("mmwr_week")
        mini_df = mini_df.groupby(["mmwr_year", "mmwr_week"]).agg(sum)
        mini_df = mini_df.reset_index()

        mini_df["bin"] = bins[0:(len(mini_df))]

        mini_df = mini_df.sort_values(by = "mmwr_week")
    

        if i == 0:
            plotting_df = mini_df
            
            i += 1

        else:
            plotting_df = pd.concat([plotting_df, mini_df])
        
    plotting_df = plotting_df.fillna(0)

    plotting_df["mmwr_year"]= pd.to_numeric(plotting_df["mmwr_year"])
    plotting_df["mmwr_week"]= pd.to_numeric(plotting_df["mmwr_week"])
    plotting_df["invasive_pneumococcal"]= pd.to_numeric(plotting_df["invasive_pneumococcal"])


    plotting_df = plotting_df.sort_values(by = ["mmwr_year", "mmwr_week"])
    plotting_x = list(range(1,len(plotting_df) + 1))

    plotting_df["x"] = plotting_x

    from sklearn.linear_model import LinearRegression

    plotting_df["COVID1"] = plotting_df["mmwr_year"] == 2019
    plotting_df["COVID2"] = plotting_df["mmwr_year"] == 2020
    plotting_df["COVID3"] = plotting_df["mmwr_year"] >= 2021


    X_data = plotting_df[["mmwr_week", "COVID1", "COVID2", "COVID3"]]
    Y_data = plotting_df["invasive_pneumococcal"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)

    binned = plotting_df.groupby(["mmwr_year","bin"]).agg(sum)
    binned = binned.reset_index()

    binned["COVID1"] = (binned["mmwr_year"] >= 2020) * 1

    X_data = binned[["bin", "COVID1"]]
    Y_data = binned["invasive_pneumococcal"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)



    modifier = {}

    for item in np.unique(bins):
        mini_binned = binned[binned["bin"] == item]
        X_data = mini_binned[["COVID1"]]
        Y_data = mini_binned["invasive_pneumococcal"]

        regressor = LinearRegression()

        regressor.fit(X_data, Y_data)
        # Now add coefficients where needed 
        modifier.update({item:regressor.coef_})

    adjustment_column = []
    for row in binned.iterrows():
        current_bin = list(row[1])[1]

        if row[1]["COVID1"] == 1:
            mod_value = list(modifier[current_bin])[0]
        
        else:
            mod_value = 0

        adjustment_column.append(mod_value)

    binned["modifier"] = adjustment_column
    binned["modifier"] = binned["modifier"] * -1
    binned["covid_adjusted"] = binned["invasive_pneumococcal"] + binned["modifier"]

    average_plot = binned.groupby("bin").agg(np.mean).reset_index()


    binned['scaled_v3'] = ((binned['covid_adjusted'] - binned.groupby('bin')['covid_adjusted'].transform('mean')) / \
                    (binned.groupby('bin')['covid_adjusted'].transform('std')))

    binned.to_csv(STATE +"_CDC_pneu.csv")

    codes = ["s5s8-d82d","4y34-2pku","xuah-ug7z","8n2k-mkiw","8rsa-pnhx","i43m-djm6","9kbf-icdi","n322-ce6f"]

    final_datasources_ip = my_dictionary()

    client = Socrata("data.cdc.gov", None)

    for code in codes:
        results = client.get(code, where = "reporting_area='" + STATE + "'")
        results_df = pd.DataFrame.from_records(results)
        final_datasources_ip.add(code, results_df)

    bins = []
    for i in range(1, 14):
        bins = bins + [i] * 4

    bins.append(13)

    i = 0
    for key in final_datasources_ip:
        columns_df = list(final_datasources_ip[key].columns)

        selection_base = ['mmwr_year', 'mmwr_week']

        if "campylobacteriosis_current_week" in columns_df:
            additional = ["campylobacteriosis_current_week"]

        elif "campylobacteriosis_current" in columns_df:
            additional = ["campylobacteriosis_current"]

        else:
            continue

        selection_final = selection_base + additional

        mini_df = final_datasources_ip[key][selection_final]
        mini_df = mini_df.rename(columns = {additional[0]:"campylobacteriosis"})
        
        mini_df["mmwr_year"]= pd.to_numeric(mini_df["mmwr_year"])
        mini_df["mmwr_week"]= pd.to_numeric(mini_df["mmwr_week"])
        mini_df["campylobacteriosis"]= pd.to_numeric(mini_df["campylobacteriosis"])
        mini_df = mini_df.fillna(0)
        mini_df = mini_df.sort_values("mmwr_week")
        mini_df = mini_df.groupby(["mmwr_year", "mmwr_week"]).agg(sum)
        mini_df = mini_df.reset_index()

        mini_df["bin"] = bins[0:(len(mini_df))]

        mini_df = mini_df.sort_values(by = "mmwr_week")

        if i == 0:
            plotting_df = mini_df
            
            i += 1

        else:
            plotting_df = pd.concat([plotting_df, mini_df])
        
    plotting_df = plotting_df.fillna(0)

    plotting_df["mmwr_year"]= pd.to_numeric(plotting_df["mmwr_year"])
    plotting_df["mmwr_week"]= pd.to_numeric(plotting_df["mmwr_week"])
    plotting_df["campylobacteriosis"]= pd.to_numeric(plotting_df["campylobacteriosis"])


    plotting_df = plotting_df.sort_values(by = ["mmwr_year", "mmwr_week"])
    plotting_x = list(range(1,len(plotting_df) + 1))

    plotting_df["x"] = plotting_x

    from sklearn.linear_model import LinearRegression

    plotting_df["COVID1"] = plotting_df["mmwr_year"] == 2019
    plotting_df["COVID2"] = plotting_df["mmwr_year"] == 2020
    plotting_df["COVID3"] = plotting_df["mmwr_year"] >= 2021


    X_data = plotting_df[["mmwr_week", "COVID1", "COVID2", "COVID3"]]
    Y_data = plotting_df["campylobacteriosis"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)

    binned = plotting_df.groupby(["mmwr_year","bin"]).agg(sum)
    binned = binned.reset_index()

    binned["COVID1"] = (binned["mmwr_year"] >= 2020) * 1

    X_data = binned[["bin", "COVID1"]]
    Y_data = binned["campylobacteriosis"]

    regressor = LinearRegression()

    regressor.fit(X_data, Y_data)

    modifier = {}

    for item in np.unique(bins):
        mini_binned = binned[binned["bin"] == item]
        X_data = mini_binned[["COVID1"]]
        Y_data = mini_binned["campylobacteriosis"]

        regressor = LinearRegression()

        regressor.fit(X_data, Y_data)
        # Now add coefficients where needed 
        modifier.update({item:regressor.coef_})

    adjustment_column = []
    for row in binned.iterrows():
        current_bin = list(row[1])[1]

        if row[1]["COVID1"] == 1:
            mod_value = list(modifier[current_bin])[0]
        
        else:
            mod_value = 0

        adjustment_column.append(mod_value)

    binned["modifier"] = adjustment_column
    binned["modifier"] = binned["modifier"] * -1
    binned["covid_adjusted"] = binned["campylobacteriosis"] + binned["modifier"]

    average_plot = binned.groupby("bin").agg(np.mean).reset_index()


    binned['scaled_v3'] = ((binned['covid_adjusted'] - binned.groupby('bin')['covid_adjusted'].transform('mean')) / \
                    (binned.groupby('bin')['covid_adjusted'].transform('std')))

    binned.to_csv(STATE +"_CDC_camp.csv")

# Time for NEON Data
final_datasources = my_dictionary()


datasources = {"DP1.20277.001":"Microbe in water (qPCR)", "DP1.10072.001":"Small Mammal (Capture Rate)", 
               "DP4.00001.001":"Summary Weather", "DP1.20138.001":"Surface Water Cell Counts", "DP1.20278.001":"Surface Water Cell Abundance (qPCR)", 
               "DP1.10092.001":"Tick Pathogen Status", "DP1.10093.001": "Tick Abundance"}

site_df = pd.read_csv("exploratory_analysis/NEON/NEON_Field_Site_Metadata_20230309.csv")


for datasource in datasources:

    url = "http://data.neonscience.org/api/v0/products/" + datasource

    response = requests.get(url)

    response_json = response.json()

    test = response_json["data"]["siteCodes"]

    if test == None:
        response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]))
    else:
        response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]["siteCodes"]))

    stripped_df = pd.DataFrame(site_df["field_site_id"])
    stripped_df.set_index("field_site_id", inplace = True)

    months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

    for i in range(2013, 2024):
        for k in months:
            a_year = str(i) + "-" + str(k)

            stripped_df[a_year] = 0

    df_tracker = my_dictionary()
    for row in response_dataframe.iterrows():
        site = list(row)[1][0]
        dates = list(row)[1][1]
        df_tracker.add(site, dates)

    for key in df_tracker:
        for item in df_tracker[key]:
            stripped_df.at[key,item] = 1

    final_datasources.add(datasource, stripped_df)

stripped_df = pd.DataFrame(site_df["field_site_id"])
stripped_df.set_index("field_site_id", inplace = True)

months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

for i in range(2013, 2024):
    for k in months:
        a_year = str(i) + "-" + str(k)

        stripped_df[a_year] = 0

temp_df = site_df[["field_site_id", "field_site_state"]]
temp_df.set_index("field_site_id", inplace = True)

field_site = temp_df["field_site_state"]

by_state = my_dictionary()

for key in datasources:
    new_dataset = final_datasources[key]
    new_dataset = new_dataset.join(field_site)
    statebased = new_dataset.groupby("field_site_state").agg(func=sum)

    by_state.add(key, statebased)

stripped_df = pd.DataFrame(site_df["field_site_state"])
stripped_df.set_index("field_site_state", inplace = True)

months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

for i in range(2013, 2024):
    for k in months:
        a_year = str(i) + "-" + str(k)

        stripped_df[a_year] = 0

for key in datasources:
    new_dataset = by_state[key]
    new_dataset.replace({2:1, 3:1, 4:1, 5:1, 6:1}, inplace = True)

    stripped_df = stripped_df + new_dataset

neon_states = list(np.unique(site_df["field_site_state"]))

site_state_neon = {}

for state in neon_states:
    field_sites = list(site_df[site_df["field_site_state"] == state]["field_site_id"])

    site_state_neon.update({state:field_sites})

import datetime
# We need to create bins that we can join

years = list(range(2012,2023))

big_dates = []
big_bins = []
for year in years: 
    datelist = pd.date_range(start = str(year) + "-01-01", end =  str(year) + "-12-31").tolist()

    str_dates = [i.strftime("%Y-%m-%d") for i in datelist]

    incr = int(len(str_dates) / 13)
    binning_list = []

    j = 1
    k = 1
    for i in range(0, len(str_dates)):
        binning_list.append(k)

        if j == incr:
            k += 1
            j = 0
            if k == 14:
                k += -1
        
        j += 1

    big_dates = big_dates + str_dates
    big_bins = big_bins + binning_list


binning_df_dates = {"date":big_dates, "bin":big_bins}
conv_df = pd.DataFrame.from_dict(binning_df_dates)
conv_df [["year", "month", "day"]]= conv_df["date"].str.split("-", expand = True)
conv_df

for id in state_ids:
    state_id = id

    # Get all dates for each site in each state
    data_item = "DP4.00001.001"
    query_list = ["wss_daily_pres", "030.01D.wss_daily_wind", "wss_daily_temp", "wss_daily_shortRad"]
    query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance']}
    all_data = {}
    all_data.update({data_item:{}})

    reference_dates = list(final_datasources[data_item].columns)

    query_dates = {}

    for site in site_state_neon[state_id]:

        query_dates.update({site:[]})

        for row in final_datasources[data_item].reset_index()[final_datasources[data_item].reset_index()["field_site_id"] == site].iterrows():
            avail = list(row[1])
            avail.pop(0)

            i = 0
            for cond in avail:
                if cond == 1:
                    query_dates[site].append(reference_dates[i])

                i += 1

    # Now create urls to query
    # Find what files are available

    for site in query_dates:
        if len(query_dates[site]) == 0:
            # Site has no entries
            continue
        else:
            for date in query_dates[site]:

                url = "http://data.neonscience.org/api/v0/data/" + data_item + "/" + site + "/" + date

                response = requests.get(url)

                response_json = response.json()

                response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]["files"]))

                files = list(response_dataframe["name"])

                for queryitem in query_list:
                    file_link = list(filter(lambda x:queryitem in x, files ))

                    if len(file_link) == 0:
                        query_list.pop(query_list.index(queryitem))
                        continue

                    elif len(file_link) != 1:
                        file_link = list(filter(lambda x:"basic" in x, file_link))

                        if len(file_link) != 1:
                            print(file_link)
                            raise ValueError("You need to make search term more specific.")
                    
                    else:
                        file_link = file_link[0]

                        query_link = url + "/" + file_link

                        imported_dataframe = pd.read_csv(query_link)

                        if queryitem not in all_data[data_item]:

                            all_data[data_item].update({queryitem:imported_dataframe})

                        else:

                            all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])

    # Bin NEON data

    summarised_data = {}
    scaled_data = {}

    for key in all_data:

        summarised_data.update({key:{}})
        scaled_data.update({key:{}})

        for table in all_data[key]:
            all_data[key][table][["year", "month", "day"]] = all_data[key][table]["date"].str.split("-", expand = True)

            all_data[key][table] = pd.merge(all_data[key][table], conv_df, left_on = ["year", "month", "day"], right_on = ["year", "month", "day"], how = "left")

            all_data[key][table] = all_data[key][table].dropna()

            all_data[key][table] = all_data[key][table][["year", "bin"] + query_tables[table]]

            summarised_data[key].update({table:all_data[key][table].groupby(["year", "bin"]).agg(np.mean).reset_index()[(["year", "bin"] + query_tables[table])]})

            scaled_data[key].update({table:summarised_data[key][table][["year", "bin"]]})

            for tab_name in query_tables[table]:
                scaled_data[key][table][tab_name] = ((summarised_data[key][table][tab_name] - summarised_data[key][table].groupby('bin')[tab_name].transform('mean')) / \
                        (summarised_data[key][table].groupby('bin')[tab_name].transform('std')))
                
    for item in scaled_data["DP4.00001.001"]:
        scaled_data["DP4.00001.001"][item].to_csv(state_id + "_" + item + ".csv")

    # Get all dates for each site in each state
    data_item = "DP1.10072.001"
    query_list = ["mam_pertrapnight"]
    query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance'],
                    "mam_pertrapnight":["trap_code","weight"]}

    all_data = {}
    all_data.update({data_item:{}})

    reference_dates = list(final_datasources[data_item].columns)

    query_dates = {}

    for site in site_state_neon[state_id]:

        query_dates.update({site:[]})

        for row in final_datasources[data_item].reset_index()[final_datasources[data_item].reset_index()["field_site_id"] == site].iterrows():
            avail = list(row[1])
            avail.pop(0)

            i = 0
            for cond in avail:
                if cond == 1:
                    query_dates[site].append(reference_dates[i])

                i += 1

    # Now create urls to query
    # Find what files are available

    for site in query_dates:
        if len(query_dates[site]) == 0:
            # Site has no entries
            continue
        else:
            for date in query_dates[site]:

                url = "http://data.neonscience.org/api/v0/data/" + data_item + "/" + site + "/" + date

                response = requests.get(url)

                response_json = response.json()

                response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]["files"]))

                files = list(response_dataframe["name"])

                for queryitem in query_list:
                    file_link = list(filter(lambda x:queryitem in x, files ))

                    if len(file_link) != 1:
                        file_link = list(filter(lambda x:"basic" in x, file_link))

                        if len(file_link) != 1:
                            raise ValueError("You need to make search term more specific.")
                        
                        else:
                            file_link = file_link[0]

                            query_link = url + "/" + file_link

                            imported_dataframe = pd.read_csv(query_link)

                            imported_dataframe[["trap_code", "trap_conv"]] = imported_dataframe["trapStatus"].str.split("-", expand = True)

                            if queryitem not in all_data[data_item]:

                                all_data[data_item].update({queryitem:imported_dataframe})

                            else:

                                all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])
                    else:
                        file_link = file_link[0]

                        query_link = url + "/" + file_link

                        imported_dataframe = pd.read_csv(query_link)

                        if queryitem not in all_data[data_item]:

                            all_data[data_item].update({queryitem:imported_dataframe([query_tables[queryitem]] + ["trap_code", "trap_conv"])})

                        else:

                            all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe([query_tables[queryitem]] + ["trap_code", "trap_conv"])])

    query_term = query_tables["mam_pertrapnight"]
    all_data['DP1.10072.001']["mam_pertrapnight"]= all_data['DP1.10072.001']["mam_pertrapnight"][["trapStatus", "collectDate","weight","trap_code", "trap_conv"]]
    all_data['DP1.10072.001']["mam_pertrapnight"]["trap_code"] = pd.to_numeric(all_data['DP1.10072.001']["mam_pertrapnight"]["trap_code"])
    all_data['DP1.10072.001']["mam_pertrapnight"]["weight"] = pd.to_numeric(all_data['DP1.10072.001']["mam_pertrapnight"]["weight"])

    all_data["DP1.10072.001"]["mam_pertrapnight"] = all_data["DP1.10072.001"]["mam_pertrapnight"][all_data["DP1.10072.001"]["mam_pertrapnight"]["trap_code"] == 5]
    all_data["DP1.10072.001"]["mam_pertrapnight"]["trap_code"] = all_data["DP1.10072.001"]["mam_pertrapnight"]["trap_code"].replace(5, 1)

    summarised_data = {}
    scaled_data = {}

    for key in all_data:

        summarised_data.update({key:{}})
        scaled_data.update({key:{}})

        for table in all_data[key]:
            all_data[key][table][["year", "month", "day"]] = all_data[key][table]["collectDate"].str.split("-", expand = True)

            all_data[key][table] = pd.merge(all_data[key][table], conv_df, left_on = ["year", "month", "day"], right_on = ["year", "month", "day"], how = "left")

            all_data[key][table] = all_data[key][table][["year", "bin"] + query_tables[table]]

            summarised_data[key].update({table:all_data[key][table].groupby(["year", "bin"]).agg(sum).reset_index()[(["year", "bin"] + query_tables[table])]})

            scaled_data[key].update({table:summarised_data[key][table][["year", "bin"]]})

            for tab_name in query_tables[table]:
                scaled_data[key][table][tab_name] = ((summarised_data[key][table][tab_name] - summarised_data[key][table].groupby('bin')[tab_name].transform('mean')) / \
                        (summarised_data[key][table].groupby('bin')[tab_name].transform('std')))
                
    scaled_data['DP1.10072.001']["mam_pertrapnight"] = scaled_data['DP1.10072.001']["mam_pertrapnight"].dropna()

    scaled_data['DP1.10072.001']["mam_pertrapnight"].to_csv(state_id + "_mammals.csv")

    # Get all dates for each site in each state
    data_item = "DP1.20138.001"
    query_list = ["amc_cellCounts"]
    query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance'],
                    "amc_cellCounts":["analysisVolume", "totalCellCount"]}
    all_data = {}
    all_data.update({data_item:{}})

    reference_dates = list(final_datasources[data_item].columns)

    query_dates = {}

    for site in site_state_neon[state_id]:

        query_dates.update({site:[]})

        for row in final_datasources[data_item].reset_index()[final_datasources[data_item].reset_index()["field_site_id"] == site].iterrows():
            avail = list(row[1])
            avail.pop(0)

            i = 0
            for cond in avail:
                if cond == 1:
                    query_dates[site].append(reference_dates[i])

                i += 1

    # Now create urls to query
    # Find what files are available
    num_sites = 0 
    for site in query_dates:
        if len(query_dates[site]) == 0:
            # Site has no entries
            continue
        else:
            num_sites += 1
            for date in query_dates[site]:

                url = "http://data.neonscience.org/api/v0/data/" + data_item + "/" + site + "/" + date

                response = requests.get(url)

                response_json = response.json()

                response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]["files"]))

                files = list(response_dataframe["name"])

                for queryitem in query_list:
                    file_link = list(filter(lambda x:queryitem in x, files ))

                    if len(file_link) > 1:
                        file_link = list(filter(lambda x:"expanded" in x, file_link))

                        if len(file_link) != 1:
                            
                            raise ValueError("You need to make search term more specific.")
                        
                        else:
                            file_link = file_link[0]

                            query_link = url + "/" + file_link

                            imported_dataframe = pd.read_csv(query_link)

                            if queryitem not in all_data[data_item]:

                                all_data[data_item].update({queryitem:imported_dataframe})

                            else:

                                all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])

                    
                    elif len(file_link) == 0:
                        continue
                        
                    else:

                        file_link = file_link[0]

                        query_link = url + "/" + file_link

                        imported_dataframe = pd.read_csv(query_link)

                        if queryitem not in all_data[data_item]:

                            all_data[data_item].update({queryitem:imported_dataframe})

                        else:

                            all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])

    # Bin NEON data

    summarised_data = {}
    scaled_data = {}

    for key in all_data:

        summarised_data.update({key:{}})
        scaled_data.update({key:{}})

        for table in all_data[key]:
            all_data[key][table]["collectDate"] = all_data[key][table]["collectDate"].str.replace("T\d\d:\d\dZ", "", regex = True)

            all_data[key][table][["year", "month", "day"]] = all_data[key][table]["collectDate"].str.split("-", expand = True)

            all_data[key][table] = pd.merge(all_data[key][table], conv_df, left_on = ["year", "month", "day"], right_on = ["year", "month", "day"], how = "left")

            all_data[key][table] = all_data[key][table][["year", "bin"] + query_tables[table]]

            all_data[key][table]['ratio'] = all_data[key][table]["totalCellCount"] / all_data[key][table]["analysisVolume"]

            query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance'],
                    "amc_cellCounts":["ratio"]}

            summarised_data[key].update({table:all_data[key][table].groupby(["year", "bin"]).agg(np.mean).reset_index()[(["year", "bin"] + query_tables[table])]})

            scaled_data[key].update({table:summarised_data[key][table][["year", "bin"]]})

            for tab_name in query_tables[table]:
                scaled_data[key][table][tab_name] = ((summarised_data[key][table][tab_name] - summarised_data[key][table].groupby('bin')[tab_name].transform('mean')) / \
                        (summarised_data[key][table].groupby('bin')[tab_name].transform('std')))

    for item in scaled_data[data_item]:
        scaled_data[data_item][item] = scaled_data[data_item][item].dropna()
        scaled_data[data_item][item].to_csv(state_id + "_" + item + ".csv")

    # Get all dates for each site in each state
    data_item = "DP1.20278.001"
    query_list = ["mga_swGroupAbundances"]
    query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance'],
                    "mga_swGroupAbundances":["dnaSampleID", "targetTaxonGroup", "meanCopyNumber"]}

    all_data = {}
    all_data.update({data_item:{}})

    reference_dates = list(final_datasources[data_item].columns)

    query_dates = {}

    for site in site_state_neon[state_id]:

        query_dates.update({site:[]})

        for row in final_datasources[data_item].reset_index()[final_datasources[data_item].reset_index()["field_site_id"] == site].iterrows():
            avail = list(row[1])
            avail.pop(0)

            i = 0
            for cond in avail:
                if cond == 1:
                    query_dates[site].append(reference_dates[i])

                i += 1

    # Now create urls to query
    # Find what files are available
    num_sites = 0 
    for site in query_dates:
        if len(query_dates[site]) == 0:
            # Site has no entries
            continue
        else:
            num_sites += 1
            for date in query_dates[site]:

                url = "http://data.neonscience.org/api/v0/data/" + data_item + "/" + site + "/" + date

                response = requests.get(url)

                response_json = response.json()

                response_dataframe = pd.DataFrame(pd.json_normalize(response_json["data"]["files"]))

                files = list(response_dataframe["name"])

                for queryitem in query_list:
                    file_link = list(filter(lambda x:queryitem in x, files ))

                    if len(file_link) > 1:
                        file_link = list(filter(lambda x:"expanded" in x, file_link))

                        if len(file_link) != 1:
                            
                            raise ValueError("You need to make search term more specific.")
                        
                        else:
                            file_link = file_link[0]

                            query_link = url + "/" + file_link

                            imported_dataframe = pd.read_csv(query_link)

                            if queryitem not in all_data[data_item]:

                                all_data[data_item].update({queryitem:imported_dataframe})

                            else:

                                all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])

                    
                    elif len(file_link) == 0:
                        continue
                        
                    else:

                        file_link = file_link[0]

                        query_link = url + "/" + file_link

                        imported_dataframe = pd.read_csv(query_link)

                        if queryitem not in all_data[data_item]:

                            all_data[data_item].update({queryitem:imported_dataframe})

                        else:

                            all_data[data_item][queryitem] = pd.concat([all_data[data_item][queryitem], imported_dataframe])

    # Bin NEON data

    summarised_data = {}
    scaled_data = {}

    for key in all_data:

        summarised_data.update({key:{}})
        scaled_data.update({key:{}})

        for table in all_data[key]:
            all_data[key][table]["collectDate"] = all_data[key][table]["collectDate"].str.replace("T\d\d:\d\dZ", "", regex = True)

            all_data[key][table][["year", "month", "day"]] = all_data[key][table]["collectDate"].str.split("-", expand = True)

            all_data[key][table] = pd.merge(all_data[key][table], conv_df, left_on = ["year", "month", "day"], right_on = ["year", "month", "day"], how = "left")

            all_data[key][table] = all_data[key][table][["year", "bin"] + query_tables[table]]

            all_data[key][table] = all_data["DP1.20278.001"]["mga_swGroupAbundances"].pivot_table(index = ["year","bin","dnaSampleID"], columns = "targetTaxonGroup", values= "meanCopyNumber").reset_index().dropna()

            all_data[key][table]["ratio"] =   all_data[key][table]["fungi"] / all_data[key][table]["bacteria and archaea"]

            summarised_data[key].update({table:all_data[key][table].groupby(["year", "bin", "dnaSampleID"]).agg(np.mean).reset_index()[(["year", "bin", "ratio"])]})

            summarised_data[key][table] = summarised_data[key][table][["year", "bin", "ratio"]]

            summarised_data[key][table] = summarised_data[key][table].groupby(["year", "bin"]).agg(np.mean).reset_index()[(["year", "bin", "ratio"])]

            query_tables = {"wss_daily_pres":["wssStaPresMean", 'wssStaPresMinimum', 'wssStaPresMaximum', 'wssStaPresVariance'], 
                    "030.01D.wss_daily_wind":['wssWindSpeedMean', 'wssWindSpeedMinimum','wssWindSpeedMaximum', 'wssWindSpeedVariance'], 
                    "wss_daily_temp":['wssTempTripleMean', 'wssTempTripleMinimum','wssTempTripleMaximum', 'wssTempTripleVariance'], 
                    "wss_daily_shortRad":['wssShortRadMean', 'wssShortRadMinimum', 'wssShortRadMaximum','wssShortRadVariance'],
                    "mga_swGroupAbundances":["ratio"]}
            
            scaled_data[key].update({table:summarised_data[key][table][["year", "bin"]]})

            for tab_name in query_tables[table]:
                scaled_data[key][table][tab_name] = ((summarised_data[key][table][tab_name] - summarised_data[key][table].groupby('bin')[tab_name].transform('mean')) / \
                        (summarised_data[key][table].groupby('bin')[tab_name].transform('std')))

    for item in scaled_data[data_item]:
        scaled_data[data_item][item] = scaled_data[data_item][item].dropna()
        scaled_data[data_item][item].to_csv(state_id + "_" + item + ".csv")