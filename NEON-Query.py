import requests
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt

class my_dictionary(dict):
 
  # __init__ function
  def __init__(self):
    self = dict()
 
  # Function to add key:value
  def add(self, key, value):
    self[key] = value

final_datasources = my_dictionary()


datasources = {"DP1.20277.001":"Microbe in water (qPCR)", "DP1.00030.001":"Methane", "DP1.00099.001":"CO2", "DP1.10072.001":"Small Mammal (Capture Rate)", 
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