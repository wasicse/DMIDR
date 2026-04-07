#!/bin/bash

# Set variables
mkdir -p ./Dispred
# file="fonPB_M.fasta"
file="Mutants_1_Res.fasta"

CONTAINER_NAME="dispredict3.03"
IMAGE_NAME="wasicse/dispredict3.0:latest"
INPUT_FASTA="/home/mkabir3/Research/27_Dispredict3_colab/Mutation/Fasta/$file"
CONTAINER_FASTA="/opt/Dispredict3.0/example/sample.fasta"
OUTPUT_CONTAINER="/opt/Dispredict3.0/output/sample_disPred.txt"
OUTPUT_HOST_1="./Dispred/${file%.fasta}.dispred"

# Step 1: Run the container (in interactive mode, detached)
echo "[1/4] Starting container..."
docker stop "$CONTAINER_NAME" || true
docker rm "$CONTAINER_NAME" || true
docker run -ti -d --name "$CONTAINER_NAME" "$IMAGE_NAME"

# Step 2: Copy the FASTA file into the container
echo "[2/4] Copying input file to container..."
docker cp "$INPUT_FASTA" "$CONTAINER_NAME":"$CONTAINER_FASTA"

# Step 3: Execute prediction inside container
echo "[3/4] Running prediction..."
docker exec -it "$CONTAINER_NAME" bash -c "
    export PATH=\"/opt/poetry/bin:\$PATH\"
    source /opt/Dispredict3.0/.venv/bin/activate
    python /opt/Dispredict3.0/script/Dispredict3.0.py -f \"$CONTAINER_FASTA\" -o \"/opt/Dispredict3.0/output/\"
"

# Step 4: Copy output back to host
echo "[4/4] Retrieving output..."
docker cp "$CONTAINER_NAME":"$OUTPUT_CONTAINER" "$OUTPUT_HOST_1"

echo "✅ Dispredict run complete."
echo "📁 Output saved to: $OUTPUT_HOST_1"
