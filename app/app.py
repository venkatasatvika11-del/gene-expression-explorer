import csv
import streamlit as st

st.title("🧬 Gene Expression Data Explorer")

st.write("Compare gene expression between Control and Treatment samples.")

# --------------------------------------------------
# Function to read CSV
# --------------------------------------------------

def read_csv(file):
    lines = file.getvalue().decode("utf-8").splitlines()
    reader = csv.DictReader(lines)

    data = list(reader)
    columns = reader.fieldnames

    return data, columns


# --------------------------------------------------
# Upload CSV file
# --------------------------------------------------

uploaded_file = st.file_uploader(
    "Upload your gene expression CSV file",
    type=["csv"]
)

if uploaded_file is not None:

    data, columns = read_csv(uploaded_file)

    st.success("CSV file loaded successfully!")

    # First column contains gene ID
    gene_column = columns[0]

    # Find Control and Treatment columns
    control_columns = []
    treatment_columns = []

    for column in columns[1:]:
        if "control" in column.lower():
            control_columns.append(column)

        if "treatment" in column.lower() or "treated" in column.lower():
            treatment_columns.append(column)

    st.write("Control columns:", control_columns)
    st.write("Treatment columns:", treatment_columns)

    # --------------------------------------------------
    # Calculate mean expression
    # --------------------------------------------------

    results = []

    for row in data:

        gene = row[gene_column]

        control_values = []

        for column in control_columns:
            control_values.append(float(row[column]))

        treatment_values = []

        for column in treatment_columns:
            treatment_values.append(float(row[column]))

        control_mean = sum(control_values) / len(control_values)
        treatment_mean = sum(treatment_values) / len(treatment_values)

        # Check whether gene is upregulated
        if treatment_mean > control_mean:
            status = "Upregulated"
        else:
            status = "Not Upregulated"

        results.append({
            "Gene": gene,
            "Control Mean": round(control_mean, 2),
            "Treatment Mean": round(treatment_mean, 2),
            "Status": status
        })

    # --------------------------------------------------
    # Display results
    # --------------------------------------------------

    st.subheader("Gene Expression Results")

    st.table(results)

else:
    st.info("Please upload a CSV file to begin.")
