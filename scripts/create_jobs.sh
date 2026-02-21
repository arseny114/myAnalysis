#!/usr/bin/env bash
#
# =============================================================================
# Скрипт для массового запуска анализа myAnalysis на существующих reconstructed файлах
# =============================================================================
#
# Назначение:
#   - Находит все файлы rec_*.root в указанной директории
#   - Для каждого файла создаёт отдельный python-скрипт конфигурации Gaudi
#   - Генерирует и подаёт на кластер отдельный job для каждого файла
#   - Выходные файлы анализа сохраняются в отдельной директории результатов
#
# Использование:
#   ./create_jobs.sh -p PROCESS_NAME -r RECO_DIR [-o ANALYSIS_ROOT] [-c CEPCSW_ROOT] [-s 0|1]
#
# Примеры:
#   ./create_jobs.sh -p E240_qqHinvi -r /cefs/higgs/.../Hinvi/E240_qqHinvi/Combined
#   ./create_jobs.sh -p E240_qqHX -r /cefs/higgs/.../HX/E240_qqHX/Reco -s 0
#
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
#          Функция вывода справки
# ──────────────────────────────────────────────────────────────────────────────
show_help() {
    cat << EOF
Скрипт для массового запуска анализа myAnalysis на кластере

Использование:
    $(basename "$0") -p PROCESS_NAME -r RECO_DIR [OPTIONS]

Обязательные параметры:
    -p, --process     Название процесса (например, E240_qqHinvi, E240_qqHX)
    -r, --reco-dir    Путь к директории с файлами реконструкции

Опциональные параметры:
    -o, --analysis-root   Корневая директория анализа
                          (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis)
    -c, --cepcsw-root     Путь к установленному CEPCSW
                          (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial)
    -s, --submit          Подавать задания на кластер (1) или только генерировать (0)
                          (по умолчанию: 1)
    -g, --group           Группа для hep_sub
                          (по умолчанию: higgs)
    -m, --memory          Требуемая память в МБ
                          (по умолчанию: 6000)
    -n, --num-files       Максимальное количество файлов для обработки
                          (по умолчанию: 0 = все файлы)
    -h, --help            Показать эту справку

Примеры:
    # Запуск для процесса H->invisible
    $(basename "$0") -p E240_qqHinvi -r /cefs/higgs/.../Hinvi/E240_qqHinvi/Combined

    # Запуск для процесса H->inclusive без подачи на кластер (отладка)
    $(basename "$0") -p E240_qqHX -r /cefs/higgs/.../HX/E240_qqHX/Reco -s 0

    # Полный вариант со всеми параметрами
    $(basename "$0") -p E240_qqHX -r /path/to/reco -o /path/to/analysis -c /path/to/cepcsw -s 1
EOF
}

# ──────────────────────────────────────────────────────────────────────────────
#          Параметры по умолчанию
# ──────────────────────────────────────────────────────────────────────────────
ANALYSIS_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis"
CEPCSW_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial"
PROCESS_NAME=""
RECO_DIR=""
SUBMIT_JOBS=1
HEP_GROUP="higgs"
MEMORY_MB=6000
MAX_FILES=0  # 0 означает обработку всех файлов

# ──────────────────────────────────────────────────────────────────────────────
#          Парсинг аргументов командной строки
# ──────────────────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--process)
            PROCESS_NAME="$2"
            shift 2
            ;;
        -r|--reco-dir)
            RECO_DIR="$2"
            shift 2
            ;;
        -o|--analysis-root)
            ANALYSIS_ROOT="$2"
            shift 2
            ;;
        -c|--cepcsw-root)
            CEPCSW_ROOT="$2"
            shift 2
            ;;
        -s|--submit)
            SUBMIT_JOBS="$2"
            shift 2
            ;;
        -g|--group)
            HEP_GROUP="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY_MB="$2"
            shift 2
            ;;
        -n|--num-files)
            MAX_FILES="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Ошибка: неизвестный параметр '$1'"
            echo "Используйте -h для просмотра справки"
            exit 1
            ;;
    esac
done

# ──────────────────────────────────────────────────────────────────────────────
#          Проверка обязательных параметров
# ──────────────────────────────────────────────────────────────────────────────
if [ -z "$PROCESS_NAME" ]; then
    echo "Ошибка: не указано название процесса (-p)"
    echo "Используйте -h для просмотра справки"
    exit 1
fi

if [ -z "$RECO_DIR" ]; then
    echo "Ошибка: не указана директория с reconstructed файлами (-r)"
    echo "Используйте -h для просмотра справки"
    exit 1
fi

# Проверка существования директории с файлами
if [ ! -d "$RECO_DIR" ]; then
    echo "Ошибка: директория не найдена: $RECO_DIR"
    exit 1
fi

# ──────────────────────────────────────────────────────────────────────────────
#          Внутренние переменные
# ──────────────────────────────────────────────────────────────────────────────
TIMESTAMP=$(date +%m%d%H%M)
SCRIPT_DIR="${ANALYSIS_ROOT}/scripts"
JOB_DIR="${ANALYSIS_ROOT}/jobs/${PROCESS_NAME}"
RES_DIR="${ANALYSIS_ROOT}/results/${PROCESS_NAME}"
LOG_DIR="${ANALYSIS_ROOT}/logs/${PROCESS_NAME}"
RECO_FILE_PATTERN="rec_${PROCESS_NAME}_*.root"

# ──────────────────────────────────────────────────────────────────────────────
#          Подготовка окружения
# ──────────────────────────────────────────────────────────────────────────────
export PATH="/cvmfs/common.ihep.ac.cn/software/hepjob/bin:$PATH"

mkdir -p "$SCRIPT_DIR" "$JOB_DIR" "$RES_DIR" "$LOG_DIR" || {
    echo "Ошибка: не удалось создать директории" >&2
    exit 1
}

# ──────────────────────────────────────────────────────────────────────────────
#          Поиск входных файлов
# ──────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "Запуск скрипта создания заданий"
echo "======================================================================"
echo
echo "Параметры запуска:"
echo "  Процесс:          ${PROCESS_NAME}"
echo "  Директория RECO:  ${RECO_DIR}"
echo "  Директория анализа: ${ANALYSIS_ROOT}"
echo "  CEPCSW_ROOT:      ${CEPCSW_ROOT}"
echo "  Подача заданий:   ${SUBMIT_JOBS}"
echo "  Группа:           ${HEP_GROUP}"
echo "  Память:           ${MEMORY_MB} MB"
echo
echo "──────────────────────────────────────────────────────────────────────"
echo "Поиск reconstructed файлов..."
echo "Путь: ${RECO_DIR}/${RECO_FILE_PATTERN}"

mapfile -t RECO_FILES < <(ls -v "${RECO_DIR}"/${RECO_FILE_PATTERN} 2>/dev/null)

# Если установлен лимит и найдено больше файлов, обрезаем массив
if [ "$MAX_FILES" -gt 0 ] && [ ${#RECO_FILES[@]} -gt "$MAX_FILES" ]; then
    echo "Ограничение количества файлов: ${MAX_FILES} из ${#RECO_FILES[@]}"
    RECO_FILES=("${RECO_FILES[@]:0:$MAX_FILES}")
fi

if [ ${#RECO_FILES[@]} -eq 0 ]; then
    echo "ОШИБКА: Не найдено ни одного файла по шаблону ${RECO_FILE_PATTERN}"
    echo "       Проверьте путь и наличие файлов!"
    exit 1
fi

printf "Найдено файлов: %d\n" "${#RECO_FILES[@]}"
echo "──────────────────────────────────────────────────────────────────────"
echo

# ──────────────────────────────────────────────────────────────────────────────
#          Основной цикл обработки каждого файла
# ──────────────────────────────────────────────────────────────────────────────
job_counter=0
for input_file in "${RECO_FILES[@]}"; do
    # Дополнительная страховка от выхода за пределы лимита внутри цикла
    if [ "$MAX_FILES" -gt 0 ] && [ "$job_counter" -ge "$MAX_FILES" ]; then
        break
    fi

    idx=$(printf "%05d" "$job_counter")
    output_file="${RES_DIR}/ana_${PROCESS_NAME}_${idx}.root"
    
    # Генерируем конфигурационный python-скрипт для этого файла
    sed -e "s|{rec_path}|${input_file}|g" \
        -e "s|{ana_path}|${output_file}|g" \
        "${SCRIPT_DIR}/temp_ana.py" > "${JOB_DIR}/ana_${idx}.py" || {
        echo "Ошибка при генерации ana_${idx}.py" >&2
        continue
    }
    
    # ── Создание и подача задания ───────────────────────────────────────────
    if [ "$SUBMIT_JOBS" -eq 1 ]; then
        sub_script="${JOB_DIR}/sub_ana_${PROCESS_NAME}_${idx}.sh"
        cat > "$sub_script" << EOF
#!/usr/bin/env bash
# Автоматически сгенерировано $(date '+%Y-%m-%d %H:%M:%S')
cd "${CEPCSW_ROOT}" || { echo "Не удалось перейти в \${CEPCSW_ROOT}"; exit 1; }
source "${CEPCSW_ROOT}/setup.sh" || { echo "Ошибка source setup.sh"; exit 1; }
echo "Запуск анализа для файла:"
echo "  Input:  ${input_file}"
echo "  Output: ${output_file}"
time ./run.sh "${JOB_DIR}/ana_${idx}.py"
echo "Завершение задания для индекса ${idx}"
EOF
        chmod +x "$sub_script"
        
        # Подача задания на кластер
        hep_sub -wt long "$sub_script" \
            -o "${LOG_DIR}/ana_${TIMESTAMP}_${PROCESS_NAME}_${idx}.out" \
            -e "${LOG_DIR}/ana_${TIMESTAMP}_${PROCESS_NAME}_${idx}.err"
        
        printf "Задание подано: %5d   %s\n" "$job_counter" "$(basename "$input_file")"
    else
        printf "Сгенерирован скрипт (без подачи): ana_%05d.py\n" "$job_counter"
    fi
    
    ((job_counter++))
done

# ──────────────────────────────────────────────────────────────────────────────
#          Итоговая информация
# ──────────────────────────────────────────────────────────────────────────────
echo
echo "======================================================================"
echo "Всего обработано файлов: ${job_counter}"
if [ "$SUBMIT_JOBS" -eq 1 ]; then
    echo "Все задания успешно поданы на кластер"
    echo "Логи будут в: ${LOG_DIR}/*.out и *.err"
else
    echo "Режим без подачи заданий (только генерация конфигов)"
fi
echo
echo "Готовые файлы анализа появятся здесь:"
echo "  ${RES_DIR}/ana_${PROCESS_NAME}_*.root"
echo "======================================================================"
