import csv
import pandas as pd
import streamlit as st

st.title("🧬 Gene Expression Data Explorer")

st.write("Compare gene expression between Control and Treatment samples.")


def read_csv(file):
    lines = file.getvalue().decode("utf-8").splitlines()
    reader = csv.DictReader(lines)
    data = list(reader)
    columns = reader.fieldnames

    return data, columns


uploaded_file = st.file_uploader(
    "Upload your gene expression CSV file",
    type=["csv"]
)


if uploaded_file is not None:

    data, columns = read_csv(uploaded_file)

    st.success("CSV file loaded successfully!")

    gene_column = columns[0]

    control_columns = []
    treatment_columns = []

    for column in columns[1:]:

        if "control" in column.lower():
            control_columns.append(column)

        if "treatment" in column.lower() or "treated" in column.lower():
            treatment_columns.append(column)

    st.subheader("Sample Groups")

    st.write("Control columns:", control_columns)
    st.write("Treatment columns:", treatment_columns)

    results = []

    for row in data:

        gene = row[gene_column]

        control_values = []

        for column in control_columns:
            value = float(row[column])
            control_values.append(value)

        treatment_values = []

        for column in treatment_columns:
            value = float(row[column])
            treatment_values.append(value)

        control_mean = sum(control_values) / len(control_values)

        treatment_mean = sum(treatment_values) / len(treatment_values)

        if treatment_mean > control_mean:
            status = "Upregulated"
        else:
            status = "Not Upregulated"

        results.append(
            {
                "Gene": gene,
                "Control Mean": round(control_mean, 2),
                "Treatment Mean": round(treatment_mean, 2),
                "Status": status
            }
        )

    st.subheader("Gene Expression Results")

    st.table(results)

    st.subheader("Control vs Treatment Expression")

    chart_data = pd.DataFrame(results)

    chart_data = chart_data.set_index("Gene")

    st.bar_chart(
        chart_data[["Control Mean", "Treatment Mean"]]
    )

else:

    st.info("Please upload a CSV file to begin.")
