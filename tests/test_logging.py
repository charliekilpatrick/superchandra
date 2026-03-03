"""
Tests for superchandra.utils.logging (get_logger, configure_logging).
"""
import logging
import pytest

from superchandra.utils.logging import (
    PACKAGE_LOG_NAME,
    get_logger,
    configure_logging,
)


class TestPackageLogName:
    def test_is_superchandra(self):
        assert PACKAGE_LOG_NAME == "superchandra"


class TestGetLogger:
    def test_returns_package_root_when_name_none(self):
        log = get_logger(None)
        assert log.name == PACKAGE_LOG_NAME
        assert isinstance(log, logging.Logger)

    def test_returns_child_logger_when_name_given(self):
        log = get_logger("foo")
        assert log.name == "superchandra.foo"

    def test_does_not_double_prefix_if_already_prefixed(self):
        log = get_logger("superchandra.bar")
        assert log.name == "superchandra.bar"

    def test_same_instance_for_same_name(self):
        a = get_logger("baz")
        b = get_logger("baz")
        assert a is b


class TestConfigureLogging:
    def test_sets_level_on_package_logger(self):
        configure_logging(level=logging.DEBUG)
        log = logging.getLogger(PACKAGE_LOG_NAME)
        assert log.level == logging.DEBUG
        # Reset so other tests are not affected
        log.setLevel(logging.NOTSET)
        for h in log.handlers[:]:
            log.removeHandler(h)

    def test_adds_stream_handler_when_none_present(self):
        log = logging.getLogger(PACKAGE_LOG_NAME)
        # Clear any handlers from previous tests
        for h in log.handlers[:]:
            log.removeHandler(h)
        configure_logging(level=logging.INFO)
        assert len(log.handlers) >= 1
        assert any(isinstance(h, logging.StreamHandler) for h in log.handlers)

    def test_does_not_duplicate_handlers_on_second_call(self):
        log = logging.getLogger(PACKAGE_LOG_NAME)
        for h in log.handlers[:]:
            log.removeHandler(h)
        configure_logging(level=logging.INFO)
        n = len(log.handlers)
        configure_logging(level=logging.INFO)
        assert len(log.handlers) == n

    def test_returns_logger(self):
        for h in logging.getLogger(PACKAGE_LOG_NAME).handlers[:]:
            logging.getLogger(PACKAGE_LOG_NAME).removeHandler(h)
        out = configure_logging(level=logging.WARNING)
        assert out is logging.getLogger(PACKAGE_LOG_NAME)


class TestLoggingOutput:
    """Check that log messages are emitted at the expected levels."""

    def test_info_warning_error_emit_when_configured(self, caplog):
        for h in logging.getLogger(PACKAGE_LOG_NAME).handlers[:]:
            logging.getLogger(PACKAGE_LOG_NAME).removeHandler(h)
        configure_logging(level=logging.DEBUG)
        log = get_logger("test_output")
        with caplog.at_level(logging.DEBUG, logger=PACKAGE_LOG_NAME):
            log.info("info message")
            log.warning("warning message")
            log.error("error message")
        assert "info message" in caplog.text
        assert "warning message" in caplog.text
        assert "error message" in caplog.text
