#!/bin/bash

# === Settings ===
CONTAINER_NAME="dispredict3.03"
IMAGE_NAME="wasicse/dispredict3.0:latest"
FASTA_DIR="./Fasta"
OUTPUT_DIR="./Dispred"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# === Main Loop ===
for i in {2..5}; do
    file="Mutants_${i}_Res.fasta"
    INPUT_FASTA="$FASTA_DIR/$file"
    CONTAINER_FASTA="/opt/Dispredict3.0/example/sample.fasta"
    OUTPUT_CONTAINER="/opt/Dispredict3.0/output/sample_disPred.txt"
    OUTPUT_HOST="$OUTPUT_DIR/${file%.fasta}.dispred"

    # Check if input file exists
    if [[ ! -f "$INPUT_FASTA" ]]; then
        echo "⚠️ Skipping: $INPUT_FASTA not found."
        continue
    fi

    echo "🔄 Running Dispredict for $file"

    # Step 1: Run container
    echo "[1/4] Starting container..."
    docker stop "$CONTAINER_NAME" &>/dev/null || true
    docker rm "$CONTAINER_NAME" &>/dev/null || true
    docker run -ti -d --name "$CONTAINER_NAME" "$IMAGE_NAME"

    # Step 2: Copy input FASTA
    echo "[2/4] Copying input file to container..."
    docker cp "$INPUT_FASTA" "$CONTAINER_NAME":"$CONTAINER_FASTA"

    # Step 3: Run prediction
    echo "[3/4] Running prediction..."
    docker exec -it "$CONTAINER_NAME" bash -c "
        export PATH=\"/opt/poetry/bin:\$PATH\"
        source /opt/Dispredict3.0/.venv/bin/activate
        python /opt/Dispredict3.0/script/Dispredict3.0.py -f \"$CONTAINER_FASTA\" -o \"/opt/Dispredict3.0/output/\"
    "

    # Step 4: Retrieve output
    echo "[4/4] Retrieving output..."
    docker cp "$CONTAINER_NAME":"$OUTPUT_CONTAINER" "$OUTPUT_HOST"

    echo "✅ Completed: $file → $OUTPUT_HOST"

    # Stop container
    docker stop "$CONTAINER_NAME" &>/dev/null
    docker rm "$CONTAINER_NAME" &>/dev/null
done

echo "🎉 All predictions complete."
