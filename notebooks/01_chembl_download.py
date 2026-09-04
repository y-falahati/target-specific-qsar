import requests
import time
import json

url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

target_id = "CHEMBL203"
limit = 1000
offset = 0

all_activities = []

while True:

    params = {
        "target_chembl_id": target_id,
        "limit": limit,
        "offset": offset
    }

    response = requests.get(url, params=params)

    print(
        f"Request: offset={offset} | "
        f"status={response.status_code}"
    )

    response.raise_for_status()

    data = response.json()

    activities = data["activities"]

    all_activities.extend(activities)

    print(f"Records received: {len(activities)}")
    print(f"Total collected: {len(all_activities)}")

    if len(activities) < limit:
        break

    offset += limit

    time.sleep(0.2)


print("\nDownload completed.")
print(f"Total activities collected: {len(all_activities)}")

output_file = "data/raw/chembl_egfr_activities.json"

with open(output_file, "w") as f:
    json.dump(all_activities, f)

print(f"Saved to: {output_file}")