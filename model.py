from sklearn.ensemble import RandomForestClassifier
import json
import pandas as pd
import os

def json_to_df(json_path: str):
    with open(json_path) as f:
        data = json.load(f)
    df = pd.json_normalize(data["Faces"])
    return df

def df_to_json(df: pd.DataFrame, output_path: str):
    data = {"Faces": df.to_dict(orient="records")}
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)

def build_feature_set(data_dir: str):
    data = []
    for file in os.listdir(data_dir):
        if not file.endswith(".json"):
            continue
        data.append(json_to_df(os.path.join(data_dir, file)))
    data = pd.concat(data)
    return data
def train():
    data = build_feature_set("/Users/tcong/mecado-onboarding/STEPs")
    X = data.drop("label", axis=1)
    y = data["label"]
    model = RandomForestClassifier()
    model.fit(X, y)
    return model
def predict(model, json_path: str):
    df = json_to_df(json_path)
    X = df.drop("label", axis=1)
    return model.predict(X)
def predict_and_save(model, json_path: str, output_path: str):
    df = json_to_df(json_path)
    X = df.drop("label", axis=1)
    df["label"] = model.predict(X)
    df_to_json(df, output_path)
if __name__ == "__main__":
    # print(json_to_df("/Users/tcong/mecado-onboarding/STEPs/2020 Aluminium Extrusion T6.json"))
    model = train()
    predict_and_save(model, "/Users/tcong/mecado-onboarding/TestSTEPs/Fusion Model.json", "/Users/tcong/mecado-onboarding/TestSTEPs/Fusion Model.json")