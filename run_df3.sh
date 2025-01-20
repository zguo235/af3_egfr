docker run -it \
    --volume $HOME/Projects/AF3/af_input:/root/af_input \
    --volume $HOME/Projects/AF3/af_output:/root/af_output \
    --volume /home/share/af_models:/root/models \
    --volume /home/share/af_db:/root/public_databases \
    --gpus '"device=0"' \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input/fold_input.json \
    --model_dir=/root/models \
    --output_dir=/root/af_output