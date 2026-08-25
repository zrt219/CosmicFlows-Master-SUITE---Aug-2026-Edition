import pytest
import time
import socket
import threading
import http.server
import socketserver
import os
from tests.cdp_client import ChromeCDPClient

_server_thread = None
_httpd = None

def ensure_http_server(port=8000):
    global _server_thread, _httpd
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', port))
    sock.close()
    if result != 0:
        class QuietHandler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, format, *args):
                pass
        
        socketserver.TCPServer.allow_reuse_address = True
        try:
            _httpd = socketserver.TCPServer(("127.0.0.1", port), QuietHandler)
            _server_thread = threading.Thread(target=_httpd.serve_forever, daemon=True)
            _server_thread.start()
            time.sleep(0.5)
        except Exception as e:
            pass

@pytest.fixture(scope="session")
def cdp_session():
    """
    Session-scoped CDP client. Launches headless Chrome and maintains connection.
    """
    ensure_http_server(8000)
    client = ChromeCDPClient(port=9222, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    # Wait for document to be ready and canvas to be mounted
    client.wait_for_condition("document.readyState === 'complete' && !!document.querySelector('#scene canvas')", timeout=10.0)
    time.sleep(1.0)
    yield client
    client.close()
    if _httpd:
        try:
            _httpd.shutdown()
        except Exception:
            pass

@pytest.fixture(scope="function")
def cdp(cdp_session):
    """
    Function-scoped fixture wrapping session CDP client.
    Clears logs before test and checks for fatal errors.
    """
    cdp_session.clear_logs() if hasattr(cdp_session, 'clear_logs') else None
    yield cdp_session
    # Assert no uncaught JS exceptions
    if cdp_session.uncaught_exceptions:
        pytest.fail(f"Uncaught JavaScript exceptions: {cdp_session.uncaught_exceptions}")
