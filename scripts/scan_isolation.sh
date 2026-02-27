#!/usr/bin/env bash
#
# =============================================================================
# Скрипт для сканирования параметров isolationDeltaR и isolationThreshold 
# в анализе myAnalysis
# =============================================================================
#
# Назначение:
#   - Запускает анализ myAnalysis с разными значениями isolationDeltaR и isolationThreshold
#   - Для каждого значения DeltaR создаётся отдельная папка
#   - Внутри папки DeltaR сохраняются файлы для всех значений Threshold
#   - Для каждого сочетания: создаёт задания, ждёт завершения, объединяет результаты
#
# Использование:
#   ./scan_isolation.sh -p PROCESS_NAME -r RECO_DIR [OPTIONS]
#
# Примеры:
#   ./scan_isolation.sh -p E240_qqHinvi -r /cefs/higgs/.../E240_qqHinvi/Reco
#   ./scan_isolation.sh -p E240_qqHinvi -r /path/to/reco -d "0.3,0.4,0.5" -t "0.1,0.5,1.0,2.0"
#
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
#          Функция вывода справки
# ──────────────────────────────────────────────────────────────────────────────
show_help() {
  cat << EOF
Скрипт для сканирования параметров isolationDeltaR и isolationThreshold в анализе myAnalysis
Использование:
$(basename "$0") -p PROCESS_NAME -r RECO_DIR [OPTIONS]

Обязательные параметры:
  -p, --process         Название процесса (например, E240_qqHinvi, E240_qqHX)
  -r, --reco-dir        Путь к директории с файлами реконструкции

Опциональные параметры:
  -d, --delta-r         Список значений isolationDeltaR через запятую
                        (по умолчанию: 0.3,0.4,0.5)
  -t, --thresholds      Список значений isolationThreshold через запятую
                        (по умолчанию: 0.1,0.2,0.5,1.0,2.0,5.0)
  -n, --num-files       Количество файлов для обработки за один запуск
                        (по умолчанию: 100)
  -o, --analysis-root   Корневая директория анализа
                        (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis)
  -c, --cepcsw-root     Путь к установленному CEPCSW
                        (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial)
  -g, --group           Группа для hep_sub
                        (по умолчанию: higgs)
  -m, --memory          Требуемая память в МБ
                        (по умолчанию: 6000)
  -k, --keep-temp       Не удалять временные файлы результатов после объединения
                        (по умолчанию: удалять)
  -h, --help            Показать эту справку

Примеры:
  # Запуск сканирования с параметрами по умолчанию
  $(basename "$0") -p E240_qqHinvi -r /cefs/higgs/.../E240_qqHinvi/Reco

  # Запуск с custom значениями DeltaR и Threshold
  $(basename "$0") -p E240_qqHinvi -r /path/to/reco -d "0.3,0.4,0.5" -t "0.1,0.5,1.0,2.0"

  # Запуск без удаления временных файлов (для отладки)
  $(basename "$0") -p E240_qqHinvi -r /path/to/reco -k
EOF
}

# ──────────────────────────────────────────────────────────────────────────────
#          Параметры по умолчанию
# ──────────────────────────────────────────────────────────────────────────────
ANALYSIS_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis"
CEPCSW_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial"
PROCESS_NAME=""
RECO_DIR=""
DELTA_R_VALUES="0.3,0.4,0.5"
THRESHOLDS="0.1,0.2,0.5,1.0,2.0,5.0"
NUM_FILES=100
HEP_GROUP="higgs"
MEMORY_MB=6000
KEEP_TEMP=0

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
    -d|--delta-r)
      DELTA_R_VALUES="$2"
      shift 2
      ;;
    -t|--thresholds)
      THRESHOLDS="$2"
      shift 2
      ;;
    -n|--num-files)
      NUM_FILES="$2"
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
    -g|--group)
      HEP_GROUP="$2"
      shift 2
      ;;
    -m|--memory)
      MEMORY_MB="$2"
      shift 2
      ;;
    -k|--keep-temp)
      KEEP_TEMP=1
      shift 1
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

if [ ! -d "$RECO_DIR" ]; then
  echo "Ошибка: директория не найдена: $RECO_DIR"
  exit 1
fi

# Проверка существования шаблона конфигурации
TEMPLATE_FILE="${ANALYSIS_ROOT}/scripts/temp_ana.py"
if [ ! -f "$TEMPLATE_FILE" ]; then
  echo "Ошибка: не найден шаблон конфигурации: $TEMPLATE_FILE"
  exit 1
fi

# ──────────────────────────────────────────────────────────────────────────────
#          Внутренние переменные
# ──────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="${ANALYSIS_ROOT}/scripts"
SCAN_ROOT_DIR="${ANALYSIS_ROOT}/scan_results/${PROCESS_NAME}"

# ──────────────────────────────────────────────────────────────────────────────
#          Подготовка окружения
# ──────────────────────────────────────────────────────────────────────────────
export PATH="/cvmfs/common.ihep.ac.cn/software/hepjob/bin:$PATH"
mkdir -p "$SCAN_ROOT_DIR" || {
  echo "Ошибка: не удалось создать директорию для результатов сканирования" >&2
  exit 1
}

# ──────────────────────────────────────────────────────────────────────────────
#          Вспомогательные функции
# ──────────────────────────────────────────────────────────────────────────────

# Функция ожидания завершения всех заданий пользователя
wait_for_jobs() {
  local process_name="$1"
  echo "Ожидание завершения заданий для процесса ${process_name}..."
  while true; do
    # Получаем ПОСЛЕДНЮЮ строку вывода hep_q -u (сводка)
    local queue_summary
    queue_summary=$(hep_q -u 2>/dev/null | tail -n 1)
    
    # Если вывод пустой или не содержит "jobs", считаем что заданий нет
    if [ -z "$queue_summary" ] || ! echo "$queue_summary" | grep -q "jobs"; then
      echo "Все задания завершены."
      return 0
    fi
    
    # Извлекаем общее число заданий из начала строки (число перед "jobs")
    local job_count
    job_count=$(echo "$queue_summary" | grep -oE '^[0-9]+' | head -n 1)
    job_count=${job_count:-0}
    
    # Если заданий нет — выходим
    if [ "$job_count" -eq 0 ]; then
      echo "Все задания завершены."
      return 0
    fi
    
    echo "  Осталось заданий: ${job_count}... Проверка через 30 секунд."
    sleep 30
  done
}

# Функция обновления значения isolationDeltaR в шаблоне
update_deltaR_in_template() {
  local delta_r="$1"
  local temp_template="${TEMPLATE_FILE}.tmp"
  # Заменяем значение после знака '=', игнорируя пробелы и комментарии после числа
  sed "s/\(myAnalysis\.isolationDeltaR[[:space:]]*=[[:space:]]*\)[0-9.]\+/\1${delta_r}/" \
    "$TEMPLATE_FILE" > "$temp_template" || return 1
  # Атомарная замена оригинала
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# Функция обновления значения isolationThreshold в шаблоне
update_threshold_in_template() {
  local threshold="$1"
  local temp_template="${TEMPLATE_FILE}.tmp"
  # Заменяем значение после знака '=', игнорируя пробелы и комментарии после числа
  sed "s/\(myAnalysis\.isolationThreshold[[:space:]]*=[[:space:]]*\)[0-9.]\+/\1${threshold}/" \
    "$TEMPLATE_FILE" > "$temp_template" || return 1
  # Атомарная замена оригинала
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# Функция восстановления исходных значений в шаблоне
restore_template_defaults() {
  local temp_template="${TEMPLATE_FILE}.tmp"
  # Восстанавливаем default значения (0.4 для DeltaR, 2.0 для Threshold)
  sed -e "s/\(myAnalysis\.isolationDeltaR[[:space:]]*=[[:space:]]*\)[0-9.]\+/\10.4/" \
      -e "s/\(myAnalysis\.isolationThreshold[[:space:]]*=[[:space:]]*\)[0-9.]\+/\12.0/" \
      "$TEMPLATE_FILE" > "$temp_template" || return 1
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# ──────────────────────────────────────────────────────────────────────────────
#          Основной цикл сканирования
# ──────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "Запуск двойного сканирования параметров изоляции"
echo "======================================================================"
echo
echo "Параметры сканирования:"
echo "  Процесс:          ${PROCESS_NAME}"
echo "  Директория RECO:  ${RECO_DIR}"
echo "  Директория анализа: ${ANALYSIS_ROOT}"
echo "  Значения DeltaR:  ${DELTA_R_VALUES}"
echo "  Значения Threshold: ${THRESHOLDS}"
echo "  Файлов за запуск: ${NUM_FILES}"
echo
echo "Результаты будут сохранены в: ${SCAN_ROOT_DIR}"
echo "======================================================================"
echo

# Преобразуем строки значений в массивы
IFS=',' read -ra DELTA_R_ARRAY <<< "$DELTA_R_VALUES"
IFS=',' read -ra THRESHOLD_ARRAY <<< "$THRESHOLDS"

total_combinations=$((${#DELTA_R_ARRAY[@]} * ${#THRESHOLD_ARRAY[@]}))
combination_counter=0
successful_combinations=0

# Внешний цикл по isolationDeltaR
deltaR_counter=0
for delta_r in "${DELTA_R_ARRAY[@]}"; do
  deltaR_counter=$((deltaR_counter + 1))
  echo ""
  echo "======================================================================"
  echo "[DeltaR ${deltaR_counter}/${#DELTA_R_ARRAY[@]}] Обработка isolationDeltaR = ${delta_r}"
  echo "======================================================================"
  echo
  
  # Создаём директорию для этого значения DeltaR
  DELTA_R_DIR="${SCAN_ROOT_DIR}/deltaR_${delta_r}"
  mkdir -p "$DELTA_R_DIR" || {
    echo "Ошибка: не удалось создать директорию ${DELTA_R_DIR}"
    continue
  }
  echo "Директория для результатов: ${DELTA_R_DIR}"
  echo
  
  # Внутренний цикл по isolationThreshold
  threshold_counter=0
  for threshold in "${THRESHOLD_ARRAY[@]}"; do
    threshold_counter=$((threshold_counter + 1))
    combination_counter=$((combination_counter + 1))
    echo "──────────────────────────────────────────────────────────────────────"
    echo "[${combination_counter}/${total_combinations}] DeltaR=${delta_r}, Threshold=${threshold}"
    echo "──────────────────────────────────────────────────────────────────────"
    echo
    
    # 1. Обновляем шаблон конфигурации с новыми значениями параметров
    echo "Шаг 1: Обновление шаблона конфигурации..."
    if ! update_deltaR_in_template "$delta_r"; then
      echo "Ошибка: не удалось обновить значение isolationDeltaR в ${TEMPLATE_FILE}"
      continue
    fi
    if ! update_threshold_in_template "$threshold"; then
      echo "Ошибка: не удалось обновить значение isolationThreshold в ${TEMPLATE_FILE}"
      # Восстанавливаем DeltaR перед продолжением
      update_deltaR_in_template "${DELTA_R_ARRAY[0]}"
      continue
    fi
    echo "  Шаблон обновлён: isolationDeltaR = ${delta_r}, isolationThreshold = ${threshold}"
    echo
    
    # 2. Запускаем создание и подачу заданий на кластер
    echo "Шаг 2: Запуск ${NUM_FILES} заданий на кластере..."
    "${SCRIPT_DIR}/create_jobs.sh" \
      -p "${PROCESS_NAME}" \
      -r "${RECO_DIR}" \
      -n "${NUM_FILES}" \
      -o "${ANALYSIS_ROOT}" \
      -c "${CEPCSW_ROOT}" \
      -g "${HEP_GROUP}" \
      -m "${MEMORY_MB}" \
      -s 1
    if [ $? -ne 0 ]; then
      echo "Ошибка: create_jobs.sh завершился с ошибкой"
      continue
    fi
    echo
    
    # 3. Ожидаем завершения всех заданий
    echo "Шаг 3: Ожидание завершения заданий..."
    wait_for_jobs "${PROCESS_NAME}"
    echo
    
    # 4. Объединяем результаты
    echo "Шаг 4: Объединение результатов..."
    # Объединяем в стандартное место
    "${SCRIPT_DIR}/merge_results.sh" -p "${PROCESS_NAME}" -o "${ANALYSIS_ROOT}" || {
      echo "Ошибка: не удалось объединить результаты для DeltaR=${delta_r}, Threshold=${threshold}"
      continue
    }
    
    # Формируем имя файла с обоими параметрами
    MERGED_BASE="merged_${PROCESS_NAME}_deltaR${delta_r}_iso${threshold}.root"
    MERGED_PATH="${DELTA_R_DIR}/${MERGED_BASE}"
    
    # Перемещаем объединённый файл в папку DeltaR
    if [ -f "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" ]; then
      mv "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" "$MERGED_PATH"
      echo "  Объединённый файл: ${MERGED_PATH}"
      successful_combinations=$((successful_combinations + 1))
    else
      echo "Ошибка: файл ${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root не найден после merge"
      continue
    fi
    echo
    
    # 5. Очищаем временные результаты (если не указано --keep-temp)
    if [ "$KEEP_TEMP" -eq 0 ]; then
      echo "Шаг 5: Очистка временных файлов..."
      RESULTS_DIR="${ANALYSIS_ROOT}/results/${PROCESS_NAME}"
      rm -f "${RESULTS_DIR}"/ana_"${PROCESS_NAME}"_*.root 2>/dev/null
      echo "  Временные файлы удалены"
    else
      echo "Шаг 5: Пропуск очистки (--keep-temp указан)"
    fi
    echo
  done
  # Конец внутреннего цикла по Threshold
  
  echo "Завершена обработка всех Threshold для DeltaR = ${delta_r}"
  echo "Файлы сохранены в: ${DELTA_R_DIR}"
  ls -lh "${DELTA_R_DIR}"/*.root 2>/dev/null || echo "  (файлы не найдены)"
  echo
done
# Конец внешнего цикла по DeltaR

# ──────────────────────────────────────────────────────────────────────────────
#          Восстановление шаблона и итоговая информация
# ──────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "Восстановление шаблона конфигурации..."
restore_template_defaults
echo "Шаблон восстановлен к значениям по умолчанию"
echo "======================================================================"
echo

echo "======================================================================"
echo "Сканирование завершено!"
echo "======================================================================"
echo
echo "Обработано комбинаций параметров: ${combination_counter} из ${total_combinations}"
echo "Успешно завершено: ${successful_combinations}"
echo
echo "Объединённые файлы:"
find "${SCAN_ROOT_DIR}" -name "merged_*.root" -type f -exec ls -lh {} \; 2>/dev/null || echo "  (файлы не найдены)"
echo "======================================================================"
