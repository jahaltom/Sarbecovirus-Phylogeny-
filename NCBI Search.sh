#!/bin/bash


conda activate seqAnalysis
# Output file
outfile="strains.tsv"
echo -e "Taxa\tIsolate\tOrganism\tHost\tDate\tLocation" > "$outfile"
# Get list of IDs
echo "🔍 Querying NCBI...Sarbecovirus "
esearch -db nucleotide -query 'txid2509511[Organism:exp] AND "complete genome"[Title] NOT "Severe acute respiratory syndrome coronavirus 2"[Organism]' | efetch -format acc  > id_list.txt
echo "📦 Downloading GenBank records and parsing..."
# Loop through IDs and extract metadata
while read acc; do
  efetch -db nucleotide -id "$acc" -format gb 2>/dev/null | \
  sed ':a;N;$!ba;s/\n[ ]\{1,\}/ /g' | \
  awk -v RS="//" -v acc="$acc" '
    BEGIN {
      isolate="NA"; organism="NA"; host="NA"; date="NA"; location="NA"
    }

    /\/isolate=/ {
      match($0, /\/isolate="([^"]+)"/, a)
      if (a[1] != "") isolate=a[1]
    }

    /\/organism=/ {
      match($0, /\/organism="([^"]+)"/, a)
      if (a[1] != "") organism=a[1]
    }

    /\/host=/ {
      match($0, /\/host="([^"]+)"/, a)
      if (a[1] != "") host=a[1]
    }

    /\/collection_date=/ {
      match($0, /\/collection_date="([^"]+)"/, a)
      if (a[1] != "") date=a[1]
    }

    /\/geo_loc_name=/ {
      match($0, /\/geo_loc_name="([^"]+)"/, a)
      if (a[1] != "") location=a[1]
    }

    END {
      print acc "\t" isolate "\t" organism "\t" host "\t" date "\t" location
    }
  '
done < id_list.txt >> "$outfile"


echo "✅ Metadata saved to: $outfile"
rm id_list.txt






python <<EOF

import pandas as pd


# Load your file
df = pd.read_csv("strains.tsv", sep="\t")

# Remove rows with NaN Date or Host
df = df[df['Date'].notna() & df['Host'].notna()]


df= df.drop_duplicates(subset=['Organism', 'Host', 'Date', 'Location'])

# Save result
df.to_csv("strains.tsv", sep="\t", index=False)

EOF






























