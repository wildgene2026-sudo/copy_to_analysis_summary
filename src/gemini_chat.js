// js/gemini_chat.js

/**
 * Isolated Gemini Chat Logic
 * This module is completely decoupled from api.js, ui.js, and app.js.
 */
function toggleGeminiSidebar() {
  const sidebar = document.getElementById('geminiSidebar');
  if (sidebar) {
    if (sidebar.classList.contains('open')) {
      sidebar.classList.remove('open');
    } else {
      sidebar.classList.add('open');
    }
  }
}

function appendGeminiMessage(role, text) {
  const history = document.getElementById('geminiChatHistory');
  if (!history) return { innerHTML: '', replaceWith: () => {} };

  const msgDiv = document.createElement('div');
  msgDiv.className = `gemini-msg ${role === 'user' ? 'gemini-msg-user' : 'gemini-msg-ai'}`;
  
  // Basic markdown handling for bold formatting
  let formattedText = text.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>');
  msgDiv.innerHTML = formattedText;
  
  history.appendChild(msgDiv);
  history.scrollTop = history.scrollHeight;
  return msgDiv;
}

// Store decoupled chat history state
let geminiChatHistoryData = [];

async function handleGeminiSubmit() {
  const inputEl = document.getElementById('geminiInput');
  if (!inputEl) return;
  
  const text = inputEl.value.trim();
  if (!text) return;

  inputEl.value = '';
  appendGeminiMessage('user', text);

  // Construct Gemini format message
  geminiChatHistoryData.push({
    role: "user",
    parts: [{text: text}]
  });

  const loadingMsg = appendGeminiMessage('ai', '<i>Thinking...</i>');

  // Attempt to grab API key from window config or local storage
  const apiKey = (typeof window.GEMINI_API_KEY !== 'undefined' ? window.GEMINI_API_KEY : null) 
              || localStorage.getItem('GEMINI_API_KEY') 
              || prompt("Please enter your Gemini API Key to continue:");
  
  if (!apiKey) {
    loadingMsg.innerHTML = "Error: GEMINI_API_KEY is not defined. Please check config.js or input when prompted.";
    return;
  }
  
  // Optionally persist the key if user typed it into prompt
  if (apiKey !== window.GEMINI_API_KEY) {
      localStorage.setItem('GEMINI_API_KEY', apiKey);
  }

  try {
    const url = `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=${apiKey}`;
    const payload = {
      contents: geminiChatHistoryData
    };

    const response = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });

    if (!response.ok) {
        if(response.status === 400 || response.status === 403) {
            localStorage.removeItem('GEMINI_API_KEY'); // Clear dead keys
        }
        throw new Error(`API Error: ${response.status} - Invalid key or request.`);
    }

    const data = await response.json();
    const replyText = data.candidates?.[0]?.content?.parts?.[0]?.text || "No response received.";
    
    loadingMsg.innerHTML = replyText.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>').replace(/\n/g, '<br>');
    
    geminiChatHistoryData.push({
      role: "model", /* Gemini uses 'model' for Assistant role */
      parts: [{text: replyText}]
    });
    
    const history = document.getElementById('geminiChatHistory');
    if(history) history.scrollTop = history.scrollHeight;

  } catch (err) {
    loadingMsg.innerHTML = `<b>Error:</b> ${err.message}`;
  }
}

document.addEventListener('DOMContentLoaded', () => {
  const inputEl = document.getElementById('geminiInput');
  if(inputEl) {
    inputEl.addEventListener('keypress', (e) => {
      if (e.key === 'Enter') {
        handleGeminiSubmit();
      }
    });
  }
});
