for sub in FERJ VACJ SEMC LEFC PIRJ CHAF
do
	sbatch --account=def-kjerbi --time=08:00:00 --job-name=job_$sub --mem=31G --ntasks-per-node=24 --nodes=1 sub_submission.sh $sub
done
