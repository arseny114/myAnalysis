#!/usr/bin/env bash
#
# =============================================================================
# Скрипт для сканирования параметров cosConeAngle и isoMaxConeEnergy 
# в анализе myAnalysis (ILC-style изоляция по энергии в конусе)
# =============================================================================
#
# Назначение:
#   - Запускает анализ myAnalysis с разными значениями cosConeAngle и isoMaxConeEnergy
#   - Для каждого значения cosConeAngle создаётся отдельная папка
#   - Внутри папки сохраняются файлы для всех значений isoMaxConeEnergy
#   - Для каждого сочетания: создаёт задания, ждёт завершения, объединяет результаты
#   - Поддерживает фоновое выполнение (работает после отключения от SSH)
#   - Ведёт логирование в указанный файл
#
# Использование:
#   ./scan_isolation.sh -p PROCESS_NAME -r RECO_DIR [OPTIONS]
#
# Примеры:
#   ./scan_isolation.sh -p E240_qqHinvi -r /cefs/higgs/.../E240_qqHinvi/Reco
#   ./scan_isolation.sh -p E240_qqHinvi -r /path/to/reco -c "0.95,0.98,0.995" -e "1.0,2.0,5.0"
#   ./scan_isolation.sh -p E240_qqHinvi -r /path/to/reco -b -l /path/to/scan.log
#
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
#          Функция вывода справки
# ──────────────────────────────────────────────────────────────────────────────
show_help() {
  cat << EOF
Скрипт для сканирования параметров изоляции (ILC-style) в анализе myAnalysis

Сканируемые параметры:
  • cosConeAngle      — косинус полу-угла конуса изоляции (cosθ ≥ cosConeAngle)
                        Пример: 0.98 ≈ угол 11.5°, 0.995 ≈ 5.7°
  • isoMaxConeEnergy  — максимальная энергия в конусе вокруг лептона (ГэВ)
                        Если энергия в конусе > порога, лептон НЕ считается изолированным

Использование:
$(basename "$0") -p PROCESS_NAME -r RECO_DIR [OPTIONS]

Обязательные параметры:
  -p, --process         Название процесса (например, E240_qqHinvi, E240_qqHX)
  -r, --reco-dir        Путь к директории с файлами реконструкции

Опциональные параметры:
  -c, --cos-cone        Список значений cosConeAngle через запятую
                        (по умолчанию: 0.95,0.98,0.995)
  -e, --cone-energy     Список значений isoMaxConeEnergy через запятую (ГэВ)
                        (по умолчанию: 1.0,2.0,5.0)
  -n, --num-files       Количество файлов для обработки за один запуск
                        (по умолчанию: 100)
  -o, --analysis-root   Корневая директория анализа
                        (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial/Analysis/myAnalysis)
  --cepcsw-root         Путь к установленному CEPCSW
                        (по умолчанию: /cefs/higgs/kositsin/CEPCSW-tutorial)
  -g, --group           Группа для hep_sub
                        (по умолчанию: higgs)
  -m, --memory          Требуемая память в МБ
                        (по умолчанию: 6000)
  -k, --keep-temp       Не удалять временные файлы результатов после объединения
                        (по умолчанию: удалять)
  -b, --background      Запустить скрипт в фоновом режиме (работает после отключения SSH)
  -l, --log-file        Путь к файлу лога (по умолчанию: scan_<process>_<timestamp>.log в ANALYSIS_ROOT)
  -h, --help            Показать эту справку

Примеры:
  # Запуск сканирования с параметрами по умолчанию (в интерактивном режиме)
  $(basename "$0") -p E240_qqHinvi -r /cefs/higgs/.../E240_qqHinvi/Reco

  # Запуск с кастомными значениями в фоновом режиме с логированием
  $(basename "$0") -p E240_qqHinvi -r /path/to/reco -c "0.97,0.99" -e "0.5,1.0,3.0" -b -l /path/to/scan.log

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
COS_CONE_VALUES="0.95,0.98,0.995"
CONE_ENERGY_VALUES="1.0,2.0,5.0"
NUM_FILES=100
HEP_GROUP="higgs"
MEMORY_MB=6000
KEEP_TEMP=0
RUN_BACKGROUND=0
LOG_FILE=""

# ──────────────────────────────────────────────────────────────────────────────
#          Обработка фонового режима и логирования
# ──────────────────────────────────────────────────────────────────────────────

# Проверяем, есть ли флаг -b в аргументах, до их обработки
for arg in "$@"; do
    if [ "$arg" = "-b" ] || [ "$arg" = "--background" ]; then
        RUN_BACKGROUND=1
        break
    fi
done

# Если запрошен фоновый режим и мы ещё не в нём, то перезапускаем скрипт через nohup
if [ "$RUN_BACKGROUND" -eq 1 ] && [ -z "${_RUNNING_IN_BACKGROUND:-}" ]; then
    # Устанавливаем путь к лог-файлу по умолчанию, если не указан
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    if [ -z "$LOG_FILE" ]; then
        LOG_FILE="${ANALYSIS_ROOT}/logs/scan_${PROCESS_NAME}_${TIMESTAMP}.log"
    fi
    
    # Создаём директорию для лога
    LOG_DIR=$(dirname "$LOG_FILE")
    mkdir -p "$LOG_DIR" 2>/dev/null
    
    echo "Запуск в фоновом режиме..."
    echo "Лог будет записываться в: $LOG_FILE"
    
    # Сохраняем PID в файл для удобства
    PID_FILE="${LOG_FILE%.log}.pid"
    
    # Собираем аргументы БЕЗ флага -b/--background и -l/--log-file, чтобы избежать рекурсии
    REEXEC_ARGS=()
    skip_next=0
    for arg in "$@"; do
        if [ "$skip_next" -eq 1 ]; then
            skip_next=0
            continue
        fi
        case "$arg" in
            -b|--background)
                # Пропускаем этот флаг
                ;;
            -l|--log-file)
                # Пропускаем флаг и его значение
                skip_next=1
                ;;
            *)
                REEXEC_ARGS+=("$arg")
                ;;
        esac
    done
    
    # Экспортируем маркер, что мы уже в фоновом режиме (на всякий случай)
    export _RUNNING_IN_BACKGROUND=1
    
    # Перезапускаем скрипт через nohup в фоне (с отфильтрованными аргументами)
    nohup "$0" "${REEXEC_ARGS[@]}" --cepcsw-root "$CEPCSW_ROOT" > "$LOG_FILE" 2>&1 &
    BG_PID=$!
    
    echo "Скрипт запущен в фоне с PID: $BG_PID"
    echo "Для просмотра лога в реальном времени: tail -f $LOG_FILE"
    echo "PID сохранён в: $PID_FILE"
    echo "$BG_PID" > "$PID_FILE"
    exit 0
fi

# Если мы уже в фоновом режиме, то помечаем это (для дочерних процессов)
if [ "$RUN_BACKGROUND" -eq 1 ]; then
    export _RUNNING_IN_BACKGROUND=1
fi

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
    -c|--cos-cone)
      COS_CONE_VALUES="$2"
      shift 2
      ;;
    -e|--cone-energy)
      CONE_ENERGY_VALUES="$2"
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
    --cepcsw-root)
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
    -b|--background)
      RUN_BACKGROUND=1
      shift 1
      ;;
    -l|--log-file)
      LOG_FILE="$2"
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

# Функция логирования с временной меткой
log_msg() {
  local level="$1"
  shift
  local msg="$*"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$timestamp] [$level] $msg"
}

# Функция ожидания завершения всех заданий пользователя
wait_for_jobs() {
  local process_name="$1"
  log_msg "INFO" "Ожидание завершения заданий для процесса ${process_name}..."
  while true; do
    local queue_summary
    queue_summary=$(hep_q -u 2>/dev/null | tail -n 1)
    
    if [ -z "$queue_summary" ] || ! echo "$queue_summary" | grep -q "jobs"; then
      log_msg "INFO" "Все задания завершены."
      return 0
    fi
    
    local job_count
    job_count=$(echo "$queue_summary" | grep -oE '^[0-9]+' | head -n 1)
    job_count=${job_count:-0}
    
    if [ "$job_count" -eq 0 ]; then
      log_msg "INFO" "Все задания завершены."
      return 0
    fi
    
    log_msg "INFO" "  Осталось заданий: ${job_count}... Проверка через 30 секунд."
    sleep 30
  done
}

# Функция обновления значения cosConeAngle в шаблоне
update_cosConeAngle_in_template() {
  local cos_val="$1"
  local temp_template="${TEMPLATE_FILE}.tmp"
  sed "s/\(myAnalysis\.cosConeAngle[[:space:]]*=[[:space:]]*\)[0-9.]\+/\1${cos_val}/" \
    "$TEMPLATE_FILE" > "$temp_template" || return 1
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# Функция обновления значения isoMaxConeEnergy в шаблоне
update_coneEnergy_in_template() {
  local energy_val="$1"
  local temp_template="${TEMPLATE_FILE}.tmp"
  sed "s/\(myAnalysis\.isoMaxConeEnergy[[:space:]]*=[[:space:]]*\)[0-9.]\+/\1${energy_val}/" \
    "$TEMPLATE_FILE" > "$temp_template" || return 1
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# Функция восстановления исходных значений в шаблоне
restore_template_defaults() {
  local temp_template="${TEMPLATE_FILE}.tmp"
  sed -e "s/\(myAnalysis\.cosConeAngle[[:space:]]*=[[:space:]]*\)[0-9.]\+/\10.98/" \
      -e "s/\(myAnalysis\.isoMaxConeEnergy[[:space:]]*=[[:space:]]*\)[0-9.]\+/\12.0/" \
      "$TEMPLATE_FILE" > "$temp_template" || return 1
  mv "$temp_template" "$TEMPLATE_FILE" || return 1
  return 0
}

# ──────────────────────────────────────────────────────────────────────────────
#          Основной цикл сканирования
# ──────────────────────────────────────────────────────────────────────────────
log_msg "INFO" "======================================================================"
log_msg "INFO" "Запуск двойного сканирования параметров изоляции (ILC-style)"
log_msg "INFO" "======================================================================"
log_msg "INFO" ""
log_msg "INFO" "Параметры сканирования:"
log_msg "INFO" "  Процесс:          ${PROCESS_NAME}"
log_msg "INFO" "  Директория RECO:  ${RECO_DIR}"
log_msg "INFO" "  Директория анализа: ${ANALYSIS_ROOT}"
log_msg "INFO" "  Значения cosConeAngle:  ${COS_CONE_VALUES}"
log_msg "INFO" "  Значения isoMaxConeEnergy: ${CONE_ENERGY_VALUES} ГэВ"
log_msg "INFO" "  Файлов за запуск: ${NUM_FILES}"
log_msg "INFO" "  Фоновый режим:    $([ "$RUN_BACKGROUND" -eq 1 ] && echo 'Да' || echo 'Нет')"
log_msg "INFO" "  Лог-файл:         ${LOG_FILE}"
log_msg "INFO" ""
log_msg "INFO" "Результаты будут сохранены в: ${SCAN_ROOT_DIR}"
log_msg "INFO" "======================================================================"
log_msg "INFO" ""

# Преобразуем строки значений в массивы
IFS=',' read -ra COS_CONE_ARRAY <<< "$COS_CONE_VALUES"
IFS=',' read -ra CONE_ENERGY_ARRAY <<< "$CONE_ENERGY_VALUES"

total_combinations=$((${#COS_CONE_ARRAY[@]} * ${#CONE_ENERGY_ARRAY[@]}))
combination_counter=0
successful_combinations=0

# Внешний цикл по cosConeAngle
cos_counter=0
for cos_val in "${COS_CONE_ARRAY[@]}"; do
  cos_counter=$((cos_counter + 1))
  angle_deg=$(python3 -c "import math; print(f'{math.degrees(math.acos(${cos_val})):.1f}')" 2>/dev/null || echo "?")
  
  log_msg "INFO" ""
  log_msg "INFO" "======================================================================"
  log_msg "INFO" "[cosConeAngle ${cos_counter}/${#COS_CONE_ARRAY[@]}] Обработка cosConeAngle = ${cos_val} (~${angle_deg}°)"
  log_msg "INFO" "======================================================================"
  log_msg "INFO" ""
  
  COS_DIR="${SCAN_ROOT_DIR}/cosCone_${cos_val}"
  mkdir -p "$COS_DIR" || {
    log_msg "ERROR" "Не удалось создать директорию ${COS_DIR}"
    continue
  }
  log_msg "INFO" "Директория для результатов: ${COS_DIR}"
  log_msg "INFO" ""
  
  # Внутренний цикл по isoMaxConeEnergy
  energy_counter=0
  for cone_energy in "${CONE_ENERGY_ARRAY[@]}"; do
    energy_counter=$((energy_counter + 1))
    combination_counter=$((combination_counter + 1))
    log_msg "INFO" "──────────────────────────────────────────────────────────────────────"
    log_msg "INFO" "[${combination_counter}/${total_combinations}] cosConeAngle=${cos_val}, isoMaxConeEnergy=${cone_energy} ГэВ"
    log_msg "INFO" "──────────────────────────────────────────────────────────────────────"
    log_msg "INFO" ""
    
    # 1. Обновляем шаблон конфигурации
    log_msg "INFO" "Шаг 1: Обновление шаблона конфигурации..."
    if ! update_cosConeAngle_in_template "$cos_val"; then
      log_msg "ERROR" "Не удалось обновить значение cosConeAngle в ${TEMPLATE_FILE}"
      continue
    fi
    if ! update_coneEnergy_in_template "$cone_energy"; then
      log_msg "ERROR" "Не удалось обновить значение isoMaxConeEnergy в ${TEMPLATE_FILE}"
      update_cosConeAngle_in_template "${COS_CONE_ARRAY[0]}"
      continue
    fi
    log_msg "INFO" "  Шаблон обновлён: cosConeAngle = ${cos_val}, isoMaxConeEnergy = ${cone_energy} ГэВ"
    log_msg "INFO" ""
    
    # 2. Запускаем создание заданий на кластер
    log_msg "INFO" "Шаг 2: Запуск ${NUM_FILES} заданий на кластере..."
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
      log_msg "ERROR" "create_jobs.sh завершился с ошибкой"
      continue
    fi
    log_msg "INFO" ""
    
    # 3. Ожидаем завершения всех заданий
    log_msg "INFO" "Шаг 3: Ожидание завершения заданий..."
    wait_for_jobs "${PROCESS_NAME}"
    log_msg "INFO" ""
    
    # 4. Объединяем результаты
    log_msg "INFO" "Шаг 4: Объединение результатов..."
    "${SCRIPT_DIR}/merge_results.sh" -p "${PROCESS_NAME}" -o "${ANALYSIS_ROOT}" || {
      log_msg "ERROR" "Не удалось объединить результаты для cosConeAngle=${cos_val}, isoMaxConeEnergy=${cone_energy}"
      continue
    }
    
    MERGED_BASE="merged_${PROCESS_NAME}_cos${cos_val}_coneE${cone_energy}GeV.root"
    MERGED_PATH="${COS_DIR}/${MERGED_BASE}"
    
    if [ -f "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" ]; then
      mv "${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root" "$MERGED_PATH"
      log_msg "INFO" "  Объединённый файл: ${MERGED_PATH}"
      successful_combinations=$((successful_combinations + 1))
    else
      log_msg "ERROR" "Файл ${ANALYSIS_ROOT}/merged_${PROCESS_NAME}.root не найден после merge"
      continue
    fi
    log_msg "INFO" ""
    
    # 5. Очищаем временные результаты
    if [ "$KEEP_TEMP" -eq 0 ]; then
      log_msg "INFO" "Шаг 5: Очистка временных файлов..."
      RESULTS_DIR="${ANALYSIS_ROOT}/results/${PROCESS_NAME}"
      rm -f "${RESULTS_DIR}"/ana_"${PROCESS_NAME}"_*.root 2>/dev/null
      log_msg "INFO" "  Временные файлы удалены"
    else
      log_msg "INFO" "Шаг 5: Пропуск очистки (--keep-temp указан)"
    fi
    log_msg "INFO" ""
  done
  
  log_msg "INFO" "Завершена обработка всех isoMaxConeEnergy для cosConeAngle = ${cos_val}"
  log_msg "INFO" "Файлы сохранены в: ${COS_DIR}"
  ls -lh "${COS_DIR}"/*.root 2>/dev/null || log_msg "WARN" "  (файлы не найдены)"
  log_msg "INFO" ""
done

# ──────────────────────────────────────────────────────────────────────────────
#          Восстановление шаблона и итоговая информация
# ──────────────────────────────────────────────────────────────────────────────
log_msg "INFO" "======================================================================"
log_msg "INFO" "Восстановление шаблона конфигурации..."
restore_template_defaults
log_msg "INFO" "Шаблон восстановлен к значениям по умолчанию"
log_msg "INFO" "======================================================================"
log_msg "INFO" ""

log_msg "INFO" "======================================================================"
log_msg "INFO" "Сканирование завершено!"
log_msg "INFO" "======================================================================"
log_msg "INFO" ""
log_msg "INFO" "Обработано комбинаций параметров: ${combination_counter} из ${total_combinations}"
log_msg "INFO" "Успешно завершено: ${successful_combinations}"
log_msg "INFO" ""
log_msg "INFO" "Объединённые файлы:"
find "${SCAN_ROOT_DIR}" -name "merged_*.root" -type f -exec ls -lh {} \; 2>/dev/null || log_msg "WARN" "  (файлы не найдены)"
log_msg "INFO" "======================================================================"
