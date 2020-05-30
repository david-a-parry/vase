vase -i ex1.vcf --af 0.4 --ac 5 --info_filters 'FS < 1' 'DB = True' | python3 ../utils/vcf2vars.py  -  > expected_outputs/test1_out.txt

vase -i ex1.vcf --max_alt 1 --keep_filters VQSRTrancheSNP99.90to99.95 | python3 ../utils/vcf2vars.py -  > expected_outputs/test2_out.txt

vase -i ex1.vcf --csq | ../utils/vcf2vars.py - > expected_outputs/test3_out.txt

vase -i ex1.vcf --impact HIGH | ../utils/vcf2vars.py - > expected_outputs/test4_out.txt

vase -i ex1.vcf  --var_type INDEL --pass  --max_alt 1 | ../utils/vcf2vars.py > expected_outputs/test5_out.txt

vase -i ex1.vcf --ped test.ped  --de_novo | ../utils/vcf2vars.py > expected_outputs/test6_out.txt

vase -i ex1.vcf --ped test.ped  --de_novo --max_alt 1  | ../utils/vcf2vars.py > expected_outputs/test7_out.txt

vase -i ex1.vcf --ped test.ped  --de_novo --het_ab 0.25 --max_alt 1  | ../utils/vcf2vars.py > expected_outputs/test8_out.txt

vase -i ex1.vcf --ped test.ped  --biallelic  --csq | ../utils/vcf2vars.py > expected_outputs/test9_out.txt

vase -i ex1.vcf --ped test.ped  --biallelic --impact HIGH MODERATE | ../utils/vcf2vars.py > expected_outputs/test10_out.txt

vase -i ex1.vcf --ped test.ped  --biallelic --impact HIGH  | ../utils/vcf2vars.py > expected_outputs/test11_out.txt
vase --ped test2.ped -i ex1.vcf --dominant --csq | ../utils/vcf2vars.py > expected_outputs/test12_out.txt
vase -i ex1.vcf --cases Sample3 Sample2 --controls Sample1 | ../utils/vcf2vars.py > expected_outputs/test13_out.txt
