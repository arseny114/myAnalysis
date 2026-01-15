#!/bin/bash

rm code_summary.txt
output_file="code_summary.txt"

> "$output_file"

# Используем find для поиска файлов и цикл for для обработки
for file in $(find . -type f -not -path './build/*' \( -name "*.cpp" -o -name "*.h" -o -name "*.py" \)); do
    echo "${file#./}:" >> "$output_file"
    echo "---" >> "$output_file"
    cat "$file" >> "$output_file"
    echo -e "\n---\n" >> "$output_file"
done

echo "Готово! Результат записан в файл: $output_file"
