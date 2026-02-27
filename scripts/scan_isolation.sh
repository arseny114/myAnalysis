#!/usr/bin/env bash
#
# =============================================================================
# Скрипт для сканирования параметра isolationThreshold в анализе myAnalysis
# =============================================================================
#
# Назначение:
#   - Запускает анализ myAnalysis с разными значениями isolationThreshold
#   - Для каждого значения: создаёт задания, ждёт завершения, объединяет результаты
#   - Сохраняет объединённые файлы с указанием значения параметра в имени
#
# Использование:
#   ./scan_isolation.sh -p PROCESS_NAME -r RECO_DIR [OPTIONS]
#
# Примеры:
#   ./scan_isolation.sh -p E240_mmHX -r /cefs/higgs/.../E240_mmHX/Reco
#   ./scan_isolation.sh -p E240_mmHX -r /path/to/reco -t "0.1,0.2,0.5,1.0,2.0"
#
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
#          Функция вывода справки
# ──────────────────────────────────────────────────────────────────────────────
show_help() {
cat << EOF
Скрипт для сканирования параметра isolationThreshold в анализе myAnalysis
Использование:
$(basename "$0") -p PROCESS_NAME -r RECO_DIR [OPTIONS]

Обязательные параметры:
  -p, --process         Название процесса (например, E240_mmHX)
  -r, --reco-dir        Путь к директории с файлами реконструкции

Опциональные параметры:
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
$(basename "$0") -p E240_mmHX -r /cefs/higgs/.../E240_mmHX/Reco

# Запуск с custom значениями порога изоляции
$(basename "$0") -p E240_mmHX -r /path/to/reco -t "0.05,0.1,0.3,0.7,1.5"

# Запуск без удаления временных файлов (для отладки)
$(basename "$0") -p E240_mmHX -r /path/to/reco -k
EOF
}

# ──────────────────────────────────────────────────────────────────────────────
#          Параметры по умолчанию
# ──────────────────────────────────────────────────────────────────────────────
ANALYSIS_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis"
CEPCSW_ROOT="/cefs/higgs/kositsin/CEPCSW-tutorial"
PROCESS_NAME=""
RECO_DIR=""
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

# ──────────────────────────────────────────────────────────────────────────────
#          Подготовка окружения
# ──────────────────────────────────────────────────────────────────────────────
export PATH="/cvmfs/common.ihep.ac.cn/software/hepjob/bin:$PATH"

# ──────────────────────────────────────────────────────────────────────────────
#          Вспомогательные функции
# ──────────────────────────────────────────────────────────────────────────────

# Функция ожидания завершения всех заданий пользователя
wait_for_jobs() {
  local process_name="$1"
  echo "Ожидание завершения заданий для процесса ${process_name}..."

  while true; do
    # Получаем ПОСЛЕДНЮЮ строку вывода hep_q -u (сводка)
    # Пример: "24 jobs; 0 completed, 0 removed, 11 idle, 13 running, 0 held, 0 suspended"
    # Важно: используем -n 1 для tail
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

# Функция обновления значения isolationThreshold в шаблоне
update_threshold_in_template() {
  local threshold="$1"
  local temp_template="${TEMPLATE_FILE}.tmp"

  # Заменяем значение после знака '=', игнорируя пробелы и комментарии после числа
  # \(...\) — захватываем часть до числа, \1 — вставляем её обратно + новое значение
  sed "s/\(myAnalysis\.isolationThreshold[[:space:]]*=[[:space:]]*\)[0-9.]\+/\1${threshold}/" \
    "$TEMPLATE_FILE" > "$temp_template" || return 1

  # Атомарная замена оригинала
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# ──────────────────────────────────────────────────────────────────────────────
#          Основной цикл сканирования
# ──────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "Запуск сканирования параметра isolationThreshold"
echo "======================================================================"
echo
echo "Параметры сканирования:"
echo "  Процесс:          ${PROCESS_NAME}"
echo "  Директория RECO:  ${RECO_DIR}"
echo "  Директория анализа: ${ANALYSIS_ROOT}"
echo "  Значения порога:  ${THRESHOLDS}"
echo "  Файлов за запуск: ${NUM_FILES}"
echo
echo "Результаты будут сохранены в: ${ANALYSIS_ROOT}"
echo "======================================================================"
echo

# Преобразуем строку значений в массив
IFS=',' read -ra THRESHOLD_ARRAY <<< "$THRESHOLDS"

scan_counter=0
for threshold in "${THRESHOLD_ARRAY[@]}"; do
  scan_counter=$((scan_counter + 1))
  echo "──────────────────────────────────────────────────────────────────────"
  echo "[${scan_counter}/${#THRESHOLD_ARRAY[@]}] Обработка threshold = ${threshold}"
  echo "──────────────────────────────────────────────────────────────────────"
  echo

  # 1. Обновляем шаблон конфигурации с новым значением порога
  echo "Шаг 1: Обновление шаблона конфигурации..."
  if ! update_threshold_in_template "$threshold"; then
    echo "Ошибка: не удалось обновить значение isolationThreshold в ${TEMPLATE_FILE}"
    continue
  fi
  echo "  Шаблон обновлён: isolationThreshold = ${threshold}"
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
  echo

  # 3. Ожидаем завершения всех заданий
  echo "Шаг 3: Ожидание завершения заданий..."
  wait_for_jobs "${PROCESS_NAME}"
  echo

  # 4. Объединяем результаты
  echo "Шаг 4: Объединение результатов..."
  MERGED_BASE="merged_${PROCESS_NAME}_iso${threshold}.root"
  MERGED_PATH="${ANALYSIS_ROOT}/${MERGED_BASE}"

  # Объединяем
  "${SCRIPT_DIR}/merge_results.sh" -p "${PROCESS_NAME}" -o "${ANALYSIS_ROOT}" || {
    echo "Ошибка: не удалось объединить результаты для threshold=${threshold}"
    continue
  }

  # Переименовываем объединённый файл
  if [ -f "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" ]; then
    mv "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" "$MERGED_PATH"
    echo "  Объединённый файл: ${MERGED_PATH}"
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

# ──────────────────────────────────────────────────────────────────────────────
#          Итоговая информация
# ──────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "Сканирование завершено!"
echo "======================================================================"
echo
echo "Обработано значений параметра: ${#THRESHOLD_ARRAY[@]}"
echo
echo "Объединённые файлы:"
ls -lh "${ANALYSIS_ROOT}"/merged_"${PROCESS_NAME}"_iso*.root 2>/dev/null || echo "  (файлы не найдены)"
echo "======================================================================"
