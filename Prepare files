import glob
import os
import pandas as pd
import numpy as np

def PROCESS_REPEATMASKER_OUTPUTS():
    # Initialize list to store all DataFrames
    all_stats = []
    
    # Step 1: Process raw RepeatMasker .fa.out files
    for file in glob.glob("*.fa.out"):
        # Standard headers for RepeatMasker output
        CORRECT_HEADERS = [
            "Score", "Pct_div", "Pct_del", "Pct_ins", "Query", "Start_in_query", 
            "End_in_query", "(left_in_scaffold)", "Orientation", "TE_ID", 
            "Class/Family", "Start_in_TE", "End_in_TE", "(left_in_TE)", "RM_ID"
        ]
        
        try:
            # Read file (skip 3 header lines, deal with whitespace separation)
            df = pd.read_csv(file, skiprows=3, header=None, sep=r'\s+')
            df.columns = CORRECT_HEADERS
            
            # Convert TE position columns to numeric, coercing errors to NaN
            df['Start_in_TE'] = pd.to_numeric(df['Start_in_TE'], errors='coerce')
            df['End_in_TE'] = pd.to_numeric(df['End_in_TE'], errors='coerce')
            df['(left_in_TE)'] = pd.to_numeric(df['(left_in_TE)'], errors='coerce')
            
            # Fix orientation for 'C' rows
            mask = df['Orientation'] == 'C'
            df.loc[mask, ['Start_in_TE', '(left_in_TE)']] = df.loc[mask, ['(left_in_TE)', 'Start_in_TE']].values
            
            # Calculate TE length (handle potential NaNs)
            df['TE_Length'] = df['End_in_TE'] - df['Start_in_TE']
            
            
            # Split Class/Family into separate columns
            split_col = df['Class/Family'].str.split('/', expand=True)
            df['Class'] = split_col[0]
            df['Family'] = split_col[1].fillna('')
            df = df.drop(columns=['Class/Family'])
            
            # Calculate statistics for each Class
            stats_df = df.groupby('Class')['TE_Length'].agg([
                ('Median_TE_Length', 'median'),
                ('Mean_TE_Length', 'mean'),
                ('Std_Dev_TE_Length', 'std'),
                ('Q1_TE_Length', lambda x: np.percentile(x, 25)),
                ('Q3_TE_Length', lambda x: np.percentile(x, 75)),
                ('IQR_TE_Length', lambda x: np.percentile(x, 75) - np.percentile(x, 25))
            ]).reset_index()
            
            # Round numeric values
            stats_df = stats_df.round(2)
            
            # Add Species column (extract from filename)
            species = os.path.basename(file).split('.')[0]
            stats_df['Species'] = species
            
            # Reorder columns
            final_df = stats_df[['Class', 'Species', 'Median_TE_Length', 'Mean_TE_Length', 
                                'Std_Dev_TE_Length', 'Q1_TE_Length', 'Q3_TE_Length', 'IQR_TE_Length']]
            
            # Append to list for consolidation
            all_stats.append(final_df)
            
            print(f"Processed: {file}")
            
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
            continue
    
    # Combine all DataFrames
    if all_stats:
        final_output = pd.concat(all_stats, ignore_index=True)
        # Save final consolidated file
        final_output.to_csv('ALL_TES_ALL_SPECIES_CONVERTED.tsv', sep='\t', index=False)
        print(f"Final output saved to ALL_TES_ALL_SPECIES_CONVERTED.tsv")
    else:
        print("No valid files were processed.")

if __name__ == "__main__":
    PROCESS_REPEATMASKER_OUTPUTS()
