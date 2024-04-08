#! /bin/sh

read -p "Please provide the analysis database: " database
database=$(echo "$database" | tr -d '[:space:]')

echo "The default gapopen and gapextend are 10 and 3, respectively."
read -p "Do you want to override them? (yes/no): " user_input
user_input=$(echo "$user_input" | tr -d '[:space:]')


if [ "$user_input" == "yes" ] || [ "$user_input" == "y" ]; then
    read -p "Please provide a custom gapopen penalty: " gapopen
    read -p "Please provide a custom gaextend penalty: " gapextend
elif [ "$user_input" == "no" ] || [ "$user_input" == "n" ]; then
    gapopen=10
    gapextend=3
else
    echo "Invalid input: $user_input"
    exit 1
fi

smdir=$(pwd)/$(ls | grep scoring_matrix)
if [ -d "$smdir" ]; then
	for file in "$smdir"/*; do
		if [ -f "$file" ]; then     
				echo
				custom_command="orthofinder -f $database --custom-matrix $file -S diamond_custom --gapopen $gapopen --gapextend $gapextend"
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

