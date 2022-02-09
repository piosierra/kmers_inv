# Workflow to find very large inversions with kmers.

## 1- Identify large inversions based on mapping of unique k-mers 
```
./kinv.sh -1 ../assemblies/2RL/DA-416_04-2RL.fa -2 ../assemblies/2RL/DA-402_03-2RL.fa -o temp2

```

## 2- Generate bed files of inversions detected between two paths on a GFA file.
```
./inv_from_paths -g ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/funestus/graphs/output_all2RL_s100000_p90_g30004000_k30/all2RL.fa.gz.a944863.f043790.61cf8a6.smooth.fix.gfa -a DA-402_03-2RL -b DA-416_04-2RL -c 2RL -p 0.8 -s 50
```