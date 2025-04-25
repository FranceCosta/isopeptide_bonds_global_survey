# Run [ML_FA](https://github.com/VivianMonzon/FAL_prediction)


After installing ML_FA:

`conda activate mlfa`

Execute with:

```bash
python3.7 /hps/software/users/agb/research/francesco/software/FAL_prediction/ML_predict.py predict \
        --fasta_seqs output/sequences.fasta \    
        --treks_dir /hps/software/users/agb/research/francesco/software/FAL_prediction/Colab/scripts/ \
        --iupred_dir /hps/software/users/agb/research/vivian/software/ \
        --analysisfolder output/analysis \
        --resultsfolder output/results \
        --jobname idp 
```