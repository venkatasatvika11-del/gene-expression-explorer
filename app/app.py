
         import csv
import pandas as pd
import streamlit as st

# --------------------------------------------------
# Title
# --------------------------------------------------

st.title("🧬 Gene Expression Data Explorer")

st.write(
    "Compare gene expression between Control and Treatment samples."
)


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


# --------------------------------------------------
# If a file is uploaded
# --------------------------------------------------

if uploaded_file is not None:

    # Read the CSV file
    data, columns = read_csv(uploaded_file)

    st.success("CSV file loaded successfully!")


    # --------------------------------------------------
    # Identify gene column
    # --------------------------------------------------

    gene_column = columns[0]


    # --------------------------------------------------
    # Find Control and Treatment columns
    # --------------------------------------------------

    control_columns = []
    treatment_columns = []

    for column in columns[1:]:

        if "control" in column.lower():
            control_columns.append(column)

        if "treatment" in column.lower() or "treated" in column.lower():
            treatment_columns.append(column)


    # --------------------------------------------------
    # Display the columns found
    # --------------------------------------------------

    st.subheader("Sample Groups")

    st.write("Control columns:", control_columns)

    st.write("Treatment columns:", treatment_columns)


    # --------------------------------------------------
    # Calculate mean expression
    # --------------------------------------------------

    results = []

    for row in data:

        # Get gene name
        gene = row[gene_column]


        # ------------------------------
        # Get Control values
        # ------------------------------

        control_values = []

        for column in control_columns:

            value = float(row[column])

            control_values.append(value)


        # ------------------------------
        # Get Treatment values
        # ------------------------------

        treatment_values = []

        for column in treatment_columns:

            value = float(row[column])

            treatment_values.append(value)


        # ------------------------------
        # Calculate means
        # ------------------------------

        control_mean = sum(control_values) / len(control_values)

        treatment_mean = sum(treatment_values) / len(treatment_values)


        # ------------------------------
        # Compare Control and Treatment
        # ------------------------------

        if treatment_mean > control_mean:

            status = "Upregulated"

        else:

            status = "Not Upregulated"


        # ------------------------------
        # Store result
        # ------------------------------

        results.append(
            {
                "Gene": gene,
                "Control Mean": round(control_mean, 2),
                "Treatment Mean": round(treatment_mean, 2),
                "Status": status
            }
        )


    # --------------------------------------------------
    # Display results
    # --------------------------------------------------

    st.subheader("Gene Expression Results")

    st.table(results)


    # --------------------------------------------------
    # Create graph
    # --------------------------------------------------

    st.subheader("Control vs Treatment Expression")

    chart_data = pd.DataFrame(results)

    chart_data = chart_data.set_index("Gene")

    st.bar_chart(
        chart_data[
            ["Control Mean", "Treatment Mean"]
        ]
    )


else:

    st.info("Please upload a CSV file to begin.")
