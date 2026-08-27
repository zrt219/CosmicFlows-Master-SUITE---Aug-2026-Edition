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

class ThreadedHTTPServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    allow_reuse_address = True
    daemon_threads = True
    def handle_error(self, request, client_address):
        pass

def ensure_http_server(port=8000):
    global _server_thread, _httpd
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', port))
    sock.close()
    if result != 0:
        class QuietHandler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, format, *args):
                pass
        
        try:
            _httpd = ThreadedHTTPServer(("127.0.0.1", port), QuietHandler)
            _server_thread = threading.Thread(target=_httpd.serve_forever, daemon=True)
            _server_thread.start()
            time.sleep(0.5)
        except Exception:
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
    Clears logs before test and ensures page health.
    """
    cdp_session.clear_logs() if hasattr(cdp_session, 'clear_logs') else None
    is_ready = False
    try:
        if cdp_session.ws and getattr(cdp_session.ws, 'connected', False):
            is_ready = cdp_session.evaluate("!!(window.cosmicflows && window.cosmicflows.timeEngine)")
    except Exception:
        is_ready = False

    if not is_ready:
        try:
            if not cdp_session.ws or not getattr(cdp_session.ws, 'connected', False) or not cdp_session.is_port_open():
                cdp_session.start()
            cdp_session.send_cdp("Page.navigate", {"url": "http://localhost:8000"})
            cdp_session.wait_for_condition("document.readyState === 'complete' && !!(window.cosmicflows && window.cosmicflows.timeEngine)", timeout=10.0)
            time.sleep(0.5)
        except Exception:
            try:
                cdp_session.close()
                cdp_session.start()
                cdp_session.send_cdp("Page.navigate", {"url": "http://localhost:8000"})
                cdp_session.wait_for_condition("document.readyState === 'complete' && !!(window.cosmicflows && window.cosmicflows.timeEngine)", timeout=10.0)
                time.sleep(0.5)
            except Exception:
                pass

    yield cdp_session
    # Assert no uncaught JS exceptions
    if cdp_session.uncaught_exceptions:
        pytest.fail(f"Uncaught JavaScript exceptions: {cdp_session.uncaught_exceptions}")

