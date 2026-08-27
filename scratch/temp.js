
window.__globalErrors = [];
window.onerror = function(msg, url, line, col, error) {
  window.__globalErrors.push({ msg: msg, url: url, line: line, col: col, stack: error ? error.stack : '' });
};
