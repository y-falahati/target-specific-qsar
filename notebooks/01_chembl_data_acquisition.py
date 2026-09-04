import requests

url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

params = {
    "target_chembl_id": "CHEMBL203",
    "limit": 5
}

response = requests.get(url, params=params)

print("Status code:", response.status_code)

data = response.json()

print("Number of activities:", len(data["activities"]))
print("Available keys:", data.keys())

print("\nFirst activity:")
print(data["activities"][0])
print("\nActivity summary:")

for activity in data["activities"]:
    print(
        activity["molecule_chembl_id"],
        "|",
        activity["standard_type"],
        "|",
        activity["standard_value"],
        activity["standard_units"],
        "| relation:",
        activity["standard_relation"],
        "| pChEMBL:",
        activity["pchembl_value"]
    )