#!/bin/bash
# logging.sh - Bash logging library with levels and colors
#
# Usage:
#   source logging.sh
#   log_info "Starting process..."
#   log_debug "Variable x=$x"
#   log_warn "Disk space low"
#   log_error "Failed to connect"
#
# Configuration (via environment variables):
#   LOG_LEVEL   - Minimum level to display: DEBUG, INFO, WARN, ERROR (default: DEBUG)
#   LOG_COLOR   - Color mode: auto, always, never (default: auto)

#==============================================================================
# Level Definitions
#==============================================================================
declare -A _LOG_LEVELS=(
    ["DEBUG"]=0
    ["INFO"]=1
    ["WARN"]=2
    ["ERROR"]=3
)

declare -A _LOG_COLORS=(
    ["DEBUG"]="\033[36m"   # Cyan
    ["INFO"]="\033[32m"    # Green
    ["WARN"]="\033[33m"    # Yellow
    ["ERROR"]="\033[31m"   # Red
    ["RESET"]="\033[0m"
)

declare -A _LOG_LABELS=(
    ["DEBUG"]="DEBUG"
    ["INFO"]="INFO "
    ["WARN"]="WARN "
    ["ERROR"]="ERROR"
)

#==============================================================================
# Configuration with Defaults
#==============================================================================
: "${LOG_LEVEL:=DEBUG}"
: "${LOG_COLOR:=auto}"

#==============================================================================
# Internal Functions
#==============================================================================

_log_get_level_num() {
    echo "${_LOG_LEVELS[$1]:-0}"
}

_log_should_log() {
    local msg_level="$1"
    local msg_num=$(_log_get_level_num "$msg_level")
    local current_num=$(_log_get_level_num "$LOG_LEVEL")
    [[ $msg_num -ge $current_num ]]
}

_log_supports_color() {
    [[ -z "${NO_COLOR:-}" ]] && [[ -t 1 ]] && [[ "${TERM:-dumb}" != "dumb" ]]
}

_log_should_use_color() {
    case "$LOG_COLOR" in
        always) return 0 ;;
        never)  return 1 ;;
        auto|*) _log_supports_color ;;
    esac
}

_log_format_message() {
    local level="$1"
    local message="$2"
    local label="${_LOG_LABELS[$level]}"

    echo "[$label] $message"
}

_log_write() {
    local level="$1"
    local message="$2"

    # Check if this level should be logged
    _log_should_log "$level" || return 0

    local formatted=$(_log_format_message "$level" "$message")

    # Write to stdout with optional color
    if _log_should_use_color; then
        local color="${_LOG_COLORS[$level]}"
        local reset="${_LOG_COLORS[RESET]}"
        echo -e "${color}${formatted}${reset}"
    else
        echo "$formatted"
    fi
}

#==============================================================================
# Public API
#==============================================================================

log_debug() {
    _log_write "DEBUG" "$*"
}

log_info() {
    _log_write "INFO" "$*"
}

log_warn() {
    _log_write "WARN" "$*"
}

log_error() {
    _log_write "ERROR" "$*" >&2
}

# Convenience function for section headers (always INFO level)
log_section() {
    log_info "$1"
}

# Check if DEBUG level is enabled
is_debug() {
    _log_should_log "DEBUG"
}

# Run a command, suppressing output unless in DEBUG mode
# Usage: run_quiet <command> [args...]
run_quiet() {
    if is_debug; then
        "$@"
    else
        "$@" > /dev/null 2>&1
    fi
}

# Validate LOG_LEVEL on source
_log_validate_level() {
    local upper_level=$(echo "$LOG_LEVEL" | tr '[:lower:]' '[:upper:]')
    if [[ -z "${_LOG_LEVELS[$upper_level]:-}" ]]; then
        echo "Warning: Invalid LOG_LEVEL '$LOG_LEVEL'. Using DEBUG." >&2
        LOG_LEVEL="DEBUG"
    else
        LOG_LEVEL="$upper_level"
    fi
    # Export so child processes inherit the log level
    export LOG_LEVEL
}

_log_validate_level
