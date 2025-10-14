// OpenAI API Configuration
// To enable the AI-powered key passage highlighting feature, 
// set your OpenAI API key here

// IMPORTANT: Uncomment and set your API key below to enable AI features
// window.OPENAI_TOKEN = 'your-openai-api-key-here';

// Example of how to set it (remove the // to activate):
// window.OPENAI_TOKEN = 'sk-proj-abcd1234...';

// Alternative: You can also set it dynamically in the browser console:
// Open DevTools > Console and type: window.OPENAI_TOKEN = 'your-key-here'

// Security Note: 
// - Only use API keys with limited permissions
// - Never commit real API keys to version control
// - The key is only used client-side for direct OpenAI API calls

console.log('OpenAI configuration loaded.');
if (typeof window !== 'undefined' && window.OPENAI_TOKEN) {
    console.log('✅ OpenAI API key detected - AI features enabled');
} else {
    console.log('⚠️  OpenAI API key not set - AI features disabled');
    console.log('💡 To enable: Set window.OPENAI_TOKEN = "your-key-here" in console or config.js');
}