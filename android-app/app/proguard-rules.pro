# Proguard Rules for CosmicFlows-4 Android APK
-keepattributes JavascriptInterface
-keepclassmembers class * {
    @android.webkit.JavascriptInterface <methods>;
}
-keep class com.zrt.cf4.AndroidBridge { *; }
-dontwarn androidx.webkit.**
