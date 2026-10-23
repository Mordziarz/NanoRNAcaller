#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="SUPPA-compliant expression cleaner")
    parser.add_argument("--expression", required=True, help="Path to IsoQuant expression TSV")
    parser.add_argument("--output", default="clean_expression.txt", help="Output clean expression file")
    args = parser.parse_args()

    df = pd.read_csv(args.expression, sep='\t', low_memory=False)

    id_col = df.columns[0]
    
    sample_cols = list(df.columns[1:])

    with open(args.output, 'w') as f:

        f.write('\t'.join(sample_cols) + '\n')
        
        for _, row in df.iterrows():
            t_id = str(row[id_col]).strip()
            if not t_id or t_id == 'nan':
                continue
                
            vals = []
            for sc in sample_cols:
                val = row[sc]
                try:
                    num = float(val) if pd.notna(val) else 0.0
                    vals.append(str(num))
                except ValueError:
                    vals.append('0.0')
                    
            f.write(t_id + '\t' + '\t'.join(vals) + '\n')

    print(f"[SUCCESS] Expression file generated strictly following SUPPA specs with {len(df)} rows and {len(sample_cols)} samples.")

if __name__ == "__main__":
    main()
