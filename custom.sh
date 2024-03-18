#! /bin/sh
# echo 
# read -p "Do you want to run the custom matrices? (yes/no): " user_input
# user_input=$(echo "$user_input" | tr -d '[:space:]')

# if [ $user_input == "yes" ] || [$user_input == "y" ]; then

smdir=$(pwd)/$(ls | grep scoring_matrix)
if [ -d "$smdir" ]; then
	for file in "$smdir"/*; do
		if [ -f "$file" ]; then     
				echo
				custom_command="orthofinder -f fungal_dataset --custom-matrix $file -S diamond_custom --gapopen 10 --gapextend 3"
				echo ">>>>> Run commond $custom_command"
				custom_result=$($custom_command)
				echo 
				echo
				echo "Command output:"
					echo "$custom_result"
		fi
	done
else
	echo "Substitution matrix directory not find!"
fi
# else
# 	echo "Operation canceled!"
# fi	
