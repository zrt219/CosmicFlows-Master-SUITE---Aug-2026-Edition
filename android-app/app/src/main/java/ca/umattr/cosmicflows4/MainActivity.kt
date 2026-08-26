package ca.umattr.cosmicflows4

import android.annotation.SuppressLint
import android.content.Context
import android.content.Intent
import android.content.res.Configuration
import android.graphics.Color
import android.net.Uri
import android.os.Build
import android.os.Bundle
import android.os.Environment
import android.os.VibrationEffect
import android.os.Vibrator
import android.os.VibratorManager
import android.util.Base64
import android.view.View
import android.view.ViewGroup
import android.view.WindowInsetsController
import android.webkit.ConsoleMessage
import android.webkit.JavascriptInterface
import android.webkit.WebChromeClient
import android.webkit.WebResourceError
import android.webkit.WebResourceRequest
import android.webkit.WebSettings
import android.webkit.WebView
import android.webkit.WebViewClient
import android.widget.Toast
import androidx.activity.ComponentActivity
import androidx.activity.OnBackPressedCallback
import androidx.core.content.FileProvider
import androidx.core.view.WindowCompat
import java.io.File
import java.io.FileOutputStream

class MainActivity : ComponentActivity() {

    private lateinit var webView: WebView

    @SuppressLint("SetJavaScriptEnabled")
    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)

        // Enable edge-to-edge layout
        WindowCompat.setDecorFitsSystemWindows(window, false)
        window.statusBarColor = Color.TRANSPARENT
        window.navigationBarColor = Color.TRANSPARENT

        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.R) {
            window.insetsController?.setSystemBarsAppearance(0, WindowInsetsController.APPEARANCE_LIGHT_STATUS_BARS)
        }

        webView = WebView(this).apply {
            layoutParams = ViewGroup.LayoutParams(
                ViewGroup.LayoutParams.MATCH_PARENT,
                ViewGroup.LayoutParams.MATCH_PARENT
            )
            setBackgroundColor(Color.BLACK)
            setLayerType(View.LAYER_TYPE_HARDWARE, null)
        }

        with(webView.settings) {
            javaScriptEnabled = true
            domStorageEnabled = true
            databaseEnabled = true
            allowFileAccess = true
            allowContentAccess = true
            allowFileAccessFromFileURLs = true
            allowUniversalAccessFromFileURLs = true
            useWideViewPort = true
            loadWithOverviewMode = true
            setSupportZoom(false)
            builtInZoomControls = false
            displayZoomControls = false
            mediaPlaybackRequiresUserGesture = false
            cacheMode = WebSettings.LOAD_DEFAULT
            mixedContentMode = WebSettings.MIXED_CONTENT_ALWAYS_ALLOW
        }

        // Attach Native Android Bridge
        webView.addJavascriptInterface(CosmicFlowsAndroidBridge(this), "AndroidBridge")

        webView.webChromeClient = object : WebChromeClient() {
            override fun onConsoleMessage(consoleMessage: ConsoleMessage?): Boolean {
                consoleMessage?.let {
                    android.util.Log.d("CosmicFlows4-JS", "[${it.messageLevel()}] ${it.message()} -- line ${it.lineNumber()}")
                }
                return true
            }
        }

        webView.webViewClient = object : WebViewClient() {
            override fun shouldOverrideUrlLoading(view: WebView?, request: WebResourceRequest?): Boolean {
                return false
            }

            override fun onReceivedError(view: WebView?, request: WebResourceRequest?, error: WebResourceError?) {
                super.onReceivedError(view, request, error)
                android.util.Log.e("CosmicFlows4-Web", "WebView Error: ${error?.description}")
            }
        }

        webView.loadUrl("file:///android_asset/index.html")
        setContentView(webView)

        onBackPressedDispatcher.addCallback(this, object : OnBackPressedCallback(true) {
            override fun handleOnBackPressed() {
                if (webView.canGoBack()) {
                    webView.goBack()
                } else {
                    isEnabled = false
                    onBackPressedDispatcher.onBackPressed()
                }
            }
        })
    }

    override fun onConfigurationChanged(newConfig: Configuration) {
        super.onConfigurationChanged(newConfig)
    }

    override fun onResume() {
        super.onResume()
        webView.onResume()
    }

    override fun onPause() {
        super.onPause()
        webView.onPause()
    }

    override fun onDestroy() {
        super.onDestroy()
        webView.destroy()
    }

    // ==========================================
    // NATIVE JAVASCRIPT INTERFACE
    // ==========================================
    inner class CosmicFlowsAndroidBridge(private val context: Context) {

        private val vibrator: Vibrator? by lazy {
            if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.S) {
                val vm = context.getSystemService(Context.VIBRATOR_MANAGER_SERVICE) as? VibratorManager
                vm?.defaultVibrator
            } else {
                @Suppress("DEPRECATION")
                context.getSystemService(Context.VIBRATOR_SERVICE) as? Vibrator
            }
        }

        @JavascriptInterface
        fun triggerHaptic(type: String) {
            try {
                if (vibrator == null || !vibrator!!.hasVibrator()) return

                if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.Q) {
                    when (type.lowercase()) {
                        "tick" -> vibrator?.vibrate(VibrationEffect.createPredefined(VibrationEffect.EFFECT_TICK))
                        "click" -> vibrator?.vibrate(VibrationEffect.createPredefined(VibrationEffect.EFFECT_CLICK))
                        "heavy" -> vibrator?.vibrate(VibrationEffect.createPredefined(VibrationEffect.EFFECT_HEAVY_CLICK))
                        "buzz" -> vibrator?.vibrate(VibrationEffect.createOneShot(80, VibrationEffect.DEFAULT_AMPLITUDE))
                        else -> vibrator?.vibrate(VibrationEffect.createPredefined(VibrationEffect.EFFECT_TICK))
                    }
                } else {
                    @Suppress("DEPRECATION")
                    when (type.lowercase()) {
                        "tick" -> vibrator?.vibrate(10)
                        "click" -> vibrator?.vibrate(25)
                        "heavy", "buzz" -> vibrator?.vibrate(60)
                        else -> vibrator?.vibrate(15)
                    }
                }
            } catch (e: Exception) {
                android.util.Log.w("CosmicFlows4-Haptics", "Haptic error: ${e.message}")
            }
        }

        @JavascriptInterface
        fun showNativeToast(message: String) {
            runOnUiThread {
                Toast.makeText(context, message, Toast.LENGTH_SHORT).show()
            }
        }

        @JavascriptInterface
        fun exportAndShare(filename: String, mimeType: String, base64Content: String) {
            runOnUiThread {
                try {
                    val cleanBase64 = if (base64Content.contains(",")) {
                        base64Content.substringAfter(",")
                    } else {
                        base64Content
                    }
                    val data = Base64.decode(cleanBase64, Base64.DEFAULT)

                    val shareDir = File(context.cacheDir, "exports")
                    if (!shareDir.exists()) shareDir.mkdirs()
                    val shareFile = File(shareDir, filename)
                    FileOutputStream(shareFile).use { it.write(data) }

                    try {
                        val downloadsDir = Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS)
                        val publicFile = File(downloadsDir, filename)
                        FileOutputStream(publicFile).use { it.write(data) }
                        Toast.makeText(context, "Saved to Downloads: $filename", Toast.LENGTH_LONG).show()
                    } catch (e: Exception) {
                        Toast.makeText(context, "Exported: $filename", Toast.LENGTH_SHORT).show()
                    }

                    val uri: Uri = FileProvider.getUriForFile(
                        context,
                        "ca.umattr.cosmicflows4.fileprovider",
                        shareFile
                    )

                    val shareIntent = Intent(Intent.ACTION_SEND).apply {
                        type = mimeType
                        putExtra(Intent.EXTRA_STREAM, uri)
                        putExtra(Intent.EXTRA_SUBJECT, "CosmicFlows-4 Publication Export: $filename")
                        putExtra(Intent.EXTRA_TEXT, "Exported from ZRT CosmicFlows-4 Cosmography Workbench (https://cf4-five.vercel.app)")
                        addFlags(Intent.FLAG_GRANT_READ_URI_PERMISSION)
                    }

                    context.startActivity(Intent.createChooser(shareIntent, "Share $filename via"))
                } catch (e: Exception) {
                    android.util.Log.e("CosmicFlows4-Export", "Export error: ${e.message}", e)
                    Toast.makeText(context, "Export failed: ${e.message}", Toast.LENGTH_SHORT).show()
                }
            }
        }
    }
}
